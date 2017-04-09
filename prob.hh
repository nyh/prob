#undef TRACE
/*
 * Given a vector of cache hit ratio (hits per request) for each of N nodes,
 * and knowing which of them is the current node, our goal is to return a
 * random vector of C different nodes (i.e., a combination of C out of N),
 * where the goals of the random distribution are:
 *
 * 1. If we send each request to the C returned nodes, the *misses per
 *    second* of all nodes will be the same. In other words, nodes with
 *    low hit ratios will be sent less work.
 *
 * 2. We know that this node is one of the N nodes. As much as possible,
 *    without breaking goal 1, we should return this node as one of the
 *    results.
 *
 * 3. We assume that *this node* got chosen uniformly randomly among the
 *    N nodes (in other words, the client chose a coordinator node, us,
 *     uniformly, and we need to choose C nodes and forward the request
 *     to them).
 */
#include <vector>
#include <utility>
#include <random>
#include <algorithm>

#include <cassert>

#include "debug.hh"

class rand_exception {};


#if 0
// This code is an incomplete attempt at giving each *combination* of K
// items a probability (calculated from the probability of each item),
// and then drawing a random combinaton. Currently, only the case of K = N-1
// is implemented, where drawing a combination is as simple as drawing one
// item which will be left *out* of the combination.
std::vector<unsigned>
randcomb2(unsigned k, const std::vector<float>& p) {
    auto n = p.size();
    if (k > p.size()) {
        throw rand_exception();
    } else if (k == 0) {
        return std::vector<unsigned>();
    }
    std::vector<unsigned> ret;
    ret.reserve(k);
    if (k == n) {
        // NOTE: In this case we do not fulfill the desired p.
        // return an error?
        std::cerr << "randcomb2: can't match p. case 1.\n";
        for (unsigned i = 0; i < n; i++) {
            ret.push_back(i);
        }
        return ret;
    } else if (k == n - 1) {
        // In this case, there are exactly k results, each result is a
        // combination x_i which misses exactly item i, so
        //
        // 1 - p(x_i) = sum_j!=i p(x_j) = k p_i
        std::vector<std::pair<unsigned,float>> np;
        np.reserve(n);
        for (unsigned i = 0; i < n; i++) {
            // TODO: If negative, throw
            float q = 1 - k*p[i];
            if (q < 0) {
                q = 0; 
            }
            else if (q > 1) {
                q = 1; 
            }
            np.emplace_back(i, q);
        }
        auto d = randone(np).first;
        std::vector<unsigned> ret;
        ret.reserve(k);
        for (unsigned i = 0; i < n; i++) {
            if (i != d) {
                ret.push_back(i);
            }
        }
        return ret;
    }
    throw rand_exception();
}
#endif

std::vector<int> randcomb(unsigned k, const std::vector<float>& p);
std::vector<float> miss_equalizing_probablities(const std::vector<float>& hit_rates);
void clip_probabilities(std::vector<float>& p, float limit);
std::vector<float> redistribute(const std::vector<float>& p, unsigned me, unsigned k);
 
 
template<typename Node>
class combination_generator {
private:
    std::vector<float> pp;
    std::vector<Node> nodes;
    int c;
public:
    combination_generator(std::vector<float>&& pp, std::vector<Node>&& nodes, int c)
        : pp(std::move(pp)), nodes(std::move(nodes)), c(c) {
        // TODO: throw if pp.size() != nodes.size() or not 1 <= c < pp.size()
    }
    std::vector<Node> get() {
        std::vector<Node> ret;
        ret.reserve(c);
        std::vector<int> r = randcomb(c, pp);
        for (int i : r) {
            ret.push_back(nodes[i]);
        }
        return ret;
    }
};




template<typename Node>
combination_generator<Node>
miss_equalizing_combination(
    const std::vector<std::pair<Node,float>>& node_hit_rate, unsigned me, int bf)
{
    auto rf = node_hit_rate.size();

    static thread_local std::default_random_engine random_engine;

    // FIXME: don't take std::pair<node,float> but separate vectors
    std::vector<float> hit_rates;
    hit_rates.reserve(rf);
    for (auto& nh : node_hit_rate) {
        hit_rates.emplace_back(nh.second);
    }
    auto p = miss_equalizing_probablities(hit_rates);
    // When we'll ask for combinations of "bf" different nodes, probabilities
    // higher than 1/bf cannot be achieved (1/bf itsef can be achieved by
    // returning this node in every returned combination). So no matter what
    // we do, we can't actually achieve the desired probabilities. Let's
    // try for the best we can
    clip_probabilities(p, 1.0f / bf);


#ifdef TRACE
    std::cout << "desired probabilities:\n";
    for (unsigned i = 0; i < p.size(); i++) {
        std::cout << node_hit_rate[i].first << ": " << p[i] <<
            ((i==me) ? " (coordinator)\n" : "\n");
    }
#endif

    // If me >= rf, this node is NOT one of the replicas, and we just need
    // to use the probabilties for these replicas, without doing the
    // redistribution to prefer the local replica.
    if (me < rf) {
        p = redistribute(p, me, bf);
    }

#ifdef TRACE
    std::cout << "returned pp: " << p << "\n";
#endif
    std::vector<Node> nodes(rf);
    for (unsigned i = 0; i < rf; i++) {
        nodes[i] = node_hit_rate[i].first;
    }
    return combination_generator<Node>(std::move(p), std::move(nodes), bf);
}

