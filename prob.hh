/*
 * Given a vector of cache hit ratio (hits per request) for each of N nodes,
 * and knowing which of them is the current node, our goal is to return a
 * random vector of K different nodes (i.e., a combination of K out of N),
 * where the goals of the random distribution are:
 *
 * 1. If we send each request to the K returned nodes, the *misses per
 *    second* of all nodes will be the same. In other words, nodes with
 *    low hit ratios will be sent less work.
 *
 * 2. We know that this node is one of the N nodes. As much as possible,
 *    without breaking goal 1, we should return this node as one of the
 *    results.
 *
 * 3. We assume that *this node* got chosen uniformly randomly among the
 *    N nodes (in other words, the client chose a coordinator node, us,
 *    uniformly, and we need to choose K nodes and forward the request
 *    to them).
 */
#include <vector>

#include "debug.hh"

class rand_exception {};

std::vector<int> randcomb(unsigned k, const std::vector<float>& p);
std::vector<float> miss_equalizing_probablities(const std::vector<float>& hit_rates);
void clip_probabilities(std::vector<float>& p, float limit);
std::vector<float> redistribute(const std::vector<float>& p, unsigned me, unsigned k);
 
 
template<typename Node>
class combination_generator {
private:
    std::vector<float> pp;
    std::vector<Node> nodes;
    int k;
public:
    combination_generator(std::vector<float>&& pp, std::vector<Node>&& nodes, int k)
        : pp(std::move(pp)), nodes(std::move(nodes)), k(k) {
        // TODO: throw if pp.size() != nodes.size() or not 1 <= k < pp.size()
    }
    std::vector<Node> get() {
        std::vector<Node> ret;
        ret.reserve(k);
        std::vector<int> r = randcomb(k, pp);
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

