#include "prob.hh"

#include <string>
#include <iostream>
#include <random>
#include <algorithm>

#include "debug.hh"

static int tests = 0, fails = 0;

static void report(bool ok, std::string msg)
{
    ++tests;
    fails += !ok;
    std::cout << (ok ? "PASS" : "FAIL") << ": " << msg << "\n";
}


// do_reorder() cycles the node_hit_rate vector to have "me" as first, because
// that's how it looks in the original caller.
std::vector<std::pair<int,float>> do_reorder(
    std::vector<std::pair<int,float>> node_hit_rate,
    int me) {
    decltype(node_hit_rate) ret;
    ret.reserve(node_hit_rate.size());
    for (unsigned i = 0; i < node_hit_rate.size(); i++) {
        auto a = node_hit_rate[(me + i) % node_hit_rate.size()];
        // add small (1%) random perturbation
        //a.second += (drand48()-0.5)*2*0.05*a.second;
        // add a consistent perturbation different for each node
        //a.second += a.second * 0.01 * i;
        ret.push_back(a);
    }
    //std::cout << ret << "\n";
    return ret;
}

void test_hit_rates(std::vector<float> hr, unsigned CL, bool reorder=false,
        unsigned iterations = 100000) {
    std::vector<std::pair<int,float>> node_hit_rate;
    int i=0;
    for (auto f : hr) {
        node_hit_rate.emplace_back(i++, f);
    }
    int N = hr.size();

    std::cout << "N=" << N << " nodes, given hit rates:\n";
    for (auto& e : node_hit_rate) {
        std::cout << "    " << e.first << " " << e.second << "\n";
    }
    std::cout << "Asking for CL=" << CL << "\n\n";

    std::random_device r;
    std::default_random_engine random_engine(r());
    std::uniform_int_distribution<> rand_node = std::uniform_int_distribution<>(0, N-1);

    // count we sent a request to each node:
    std::vector<int> count;
    count.resize(N);
    // count of times that a particular coordinator sent to a particular node:
    std::vector<std::vector<int>> count2;
    count2.resize(N);
    for (auto& v : count2) {
        v.resize(N);
    }
    // count of *invalid* combinations drawn. Should be 0 of course
    int count_invalid = 0;

    // Do "iterations" loop iterations, in each we pick a coordinator randomly, and
    // then see which CL nodes it will choose.
    for (unsigned i = 0; i < iterations; i++) {
        int coord = rand_node(random_engine);
        std::vector<int> c;
        if (reorder) {
            auto gen = miss_equalizing_combination(do_reorder(node_hit_rate, coord), 0, CL);
            c = gen.get();
        } else {
            auto gen = miss_equalizing_combination(node_hit_rate, coord, CL);
            c = gen.get();
        }
        // sort for nicer debugging printout (more useful with uniq) but also to
        // make it easier to find invalid combinations with duplicate nodes
        std::sort(c.begin(), c.end());
        bool invalid = false;
        int prev = -1;
        for(auto& s : c) {
            if (prev == s) {
                // found duplicate in returned combination! It's easy to
                // check because we sorted the combination.
                invalid = true;
            }
            prev = s;
            count[s]++;
            count2[coord][s]++;
        }
        if (invalid) {
            std::cout << "invalid combination: ";
            for (auto& s : c) {
                std::cout << s << " ";
            }
            std::cout << "\n";
            count_invalid++;
        }
//        std::cout << "\n";
    }

    float sum = 0;
    for(int i = 0; i < N; i++) {
        sum += count[i];
    }
    std::cout << "work sent to nodes: ";
    for(int i = 0; i < N; i++) {
        std::cout << count[i]/sum << " ";
    }
    std::cout << "\n";
    std::cout << "          expected: ";
    std::vector<float> p;
    float p_sum = 0;
    for(int i = 0; i < N; i++) {
        float r = 1/(1-node_hit_rate[i].second);
        p.push_back(r);
        p_sum += r;
    }
    for(int i = 0; i < N; i++) {
        p[i] /= p_sum;
    }
    for(int i = 0; i < N; i++) {
        std::cout << p[i] << " ";
    }
    std::cout << "\n";
    std::cout << "             error: ";
    for(int i = 0; i < N; i++) {
        std::cout << (count[i]/sum/p[i] - 1)*100 << "% ";
    }
    std::cout << "\n";
    for(int i = 0; i < N; i++) {
        report(std::abs((count[i]/sum/p[i] - 1)*100) < 4, "error less than 4%");
    }
    std::cout << "misses (work * (1-hitrate)), normalized. Should be all 1:\n\t";
    float min = 1.0;
    for(int i = 0; i < N; i++) {
        min = std::min(min, (1-node_hit_rate[i].second)*count[i]/sum);
    }
    for(int i = 0; i < N; i++) {
        std::cout << (1-node_hit_rate[i].second)*count[i]/sum/min << " ";
    }
    std::cout << "\n";
    std::cout << "\n";
    for (int i = 0; i < N; i++) {
        int total = 0;
        for (int j = 0; j < N; j++) {
            total += count2[i][j];
        }
        std::cout << "Work sent by coordinator " << i << ": (of total " << total << ")\n";
        for (int j = 0; j < N; j++) {
            std:: cout << "\t" << j << " " << (float)count2[i][j]/total << "\n";
        }
        // Node should send the maximum possible to itself: NP or 1/C, whichever is
        // lower (since it cannot keep more than 1/C locally).
        float e = std::min(N * p[i], 1.0f / CL);
        report(std::abs(((float)count2[i][i]/total) - e) < 1e-2,
                "node should keep as much as possible locally");
        // Node cannot send to any node (including itself) more than 1/C
        for (int j = 0; j < N; j++) {
            report(1.0f / CL - ((float)count2[i][i]/total) > -1e-3,
                "cannot send more than 1/CL to individual node");
        }
    }
    std::cout << "\n";
    if (count_invalid) {
        std::cout << "ERROR: Found invalid combinations: " << count_invalid << "\n";
    }
    report(count_invalid == 0, "no random combination can have the same node more than once");
    std::cout << "----------------------------------------------------------\n";
}

int
main() {
    // Test starting with various hit rates (the first step with those is to calculate
    // expected probabilities, and continue from them).

    // Test for CL=1
    test_hit_rates({0.8, 0.65, 0.55}, 1);
    test_hit_rates({0.8, 0.8, 0.8}, 1);
    test_hit_rates({0.9, 0.25, 0.15}, 1);
    
    // Tests with N=3, CL=2
    test_hit_rates({0.8, 0.8, 0.8}, 2);
    test_hit_rates({0.8, 0.8, 0.2}, 2);
    test_hit_rates({0.8, 0.57, 0.8}, 2);
    test_hit_rates({0.81, 0.79, 0.57}, 2);
    test_hit_rates({0.81, 0.79, 0.27}, 2);
    test_hit_rates({0.87, 0.83, 0.75}, 2);
    test_hit_rates({0.95, 0.95, 0.15}, 2);
    test_hit_rates({0.8, 0.65, 0.55}, 2);
    // Reorder the nodes on different coordinators (cyclic reordering so the
    // coordinator is first) and see that it doesn't mess up the decisions like
    // happened in previous versions.
    test_hit_rates({0.8, 0.57, 0.8}, 2, true);

    // Tests with N=4, CL=3
    test_hit_rates({0.90, 0.89, 0.89, 0.40}, 3);
#if 0
    // BUG! This test has original probabilities 0.356234, 0.409669, 0.117048,
    // 0.117048, clipped at 1/CL to 0.333333, 0.333333, 0.166667, 0.166667
    // However, even those are not achieved, and we achieve the very wrong
    // 0.291834 0.291737 0.166724 0.249706. Note how the 3rd and 4th nodes
    // which, no matter what, should have received the same amount of work,
    // did not. So it's definitely a bug.
    test_hit_rates({0.77, 0.80, 0.30, 0.30}, 3);
#endif
#if 0
    // BUG! This test cannot fully succeed because P(0.8) = 0.40 and it should
    // be clipped to 1/CL=0.33 (and it is). However, we're getting a lot of
    // invalid combination printouts - that should NOT happen. It happens
    // because the algorithm generates in pp 0.441 which is > 1/CL, and that's
    // A BUG.
    test_hit_rates({0.77, 0.80, 0.30, 0.32}, 3);
#endif
#if 0
    // BUG! The error level is very low, but still we get a lot of invalid
    // combinations.
    // This happens because we get in pp a probability higher than 0.33 -
    // the first stage divided p node 1's causing me (which has the highest
    // deficit) to partitipate in the division twice and give it pp[1] = 0.12
    // and then later at the end we 2 remaining mixed nodes with surplus
    // 0.316, and add that to the 0.12 and get over 1/CL = 0.333....
    // Perhaps like we fixed the one-mixed-node-remaining case we also need
    // to fix this case? But it will be very messy to track for the different
    // me without running the full algorithm :-(
    test_hit_rates({0.66, 0.66, 0.34, 0.32}, 3);
#endif

    // Tests with N=5, CL=2
    test_hit_rates({0.79, 0.78, 0.77, 0.80, 0.32}, 2);

    // Tests with N=5, CL=3
    test_hit_rates({0.79, 0.78, 0.77, 0.80, 0.32}, 3);

    // Tests with N=7, CL=4
    test_hit_rates({0.79, 0.78, 0.77, 0.80, 0.33, 0.33, 0.3}, 4);

#if 0
    // In this case, it is impossible to reproduce the desired probabilities because
    // node 0 wants probability 0.652 > 0.5. So it is rather pointless to test this
    // call. However, we should separately test the clipping algorithm.
    test_hit_rates({0.95, 0.85, 0.75}, 2);
    // Here, as we saw before, one of the probabilities is 0.348 > 0.333, so we
    // can't achieve the exact probabilities.
    test_hit_rates({0.90, 0.89, 0.91, 0.40}, 3);
#endif

    std::cout << "SUMMARY: " << tests << " tests, " << fails << " failures\n";
    return fails == 0 ? 0 : 1;
}
