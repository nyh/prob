#include "prob.hh"

#include <string>

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

int
main() {
    bool reorder = false;
#if 0
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.8},
        {1, 0.8},
        {2, 0.8},
    };
    int CL = 2;
#elif 0
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.8},
        {1, 0.8},
        {2, 0.2},
    };
    int CL = 2;
#elif 0
    // Reorder the nodes on different coordinators (cyclic reordering
    // so the coordinator is first) and see that it doesn't mess up the
    // decisions like happened in previous versions.
    reorder = true;
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.8},
        {1, 0.57},
        {2, 0.8},
    };
    int CL = 2;
#elif 0
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.81},
        {1, 0.79},
        {2, 0.57},
    };
    int CL = 2;
#elif 0
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.81},
        {1, 0.79},
        {2, 0.27},
    };
    int CL = 2;
#elif 0
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.87},
        {1, 0.83},
        {2, 0.75},
    };
    int CL = 2;
#elif 0
    // This case is impossible to reproduce the desired
    // probabilities because node 0 wants probability 0.652 > 0.5.
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.95},
        {1, 0.85},
        {2, 0.75},
    };
    int CL = 2;
#elif 0
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.95},
        {1, 0.95},
        {2, 0.15},
    };
    int CL = 2;
#elif 0
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.79},
        {1, 0.78},
        {2, 0.77},
        {3, 0.80},
        {4, 0.32},
    };
    int CL = 2;
#elif 0
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.79},
        {1, 0.78},
        {2, 0.77},
        {3, 0.80},
        {4, 0.32},
    };
    int CL = 3;
#elif 0
    // Test for the algorithm fix for the case where after the algorithm we
    // we are left with just one mixed node, and need to redistribute some of
    // the pp changes.
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.80},
        {1, 0.65},
        {2, 0.55},
    };
    int CL = 2;
#elif 0
    // Test for CL=1
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.80},
        {1, 0.65},
        {2, 0.55},
    };
    int CL = 1;
#elif 0
    // Test for CL=1
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.80},
        {1, 0.80},
        {2, 0.80},
    };
    int CL = 1;
#elif 0
    // This test is not working correctly. probably because P > 1/CL?
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.79},
        {1, 0.78},
        {2, 0.77},
        {3, 0.80},
        {4, 0.33},
        {5, 0.33},
        {6, 0.30},
    };
    int CL = 4;
#elif 0
    // This test is cannot work correctly because P of 0.80 is 0.409 > 1/CL?
    // However, I would expect the probabilities to get clamped to 0.33, and
    // I see them clamped to 0.29 - why??
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.77},
        {1, 0.80},
        {2, 0.30},
        {3, 0.30},
    };
    int CL = 3;
#elif 0
    // BUG?
    // This test cannot fully succeed because P(0.8) = 0.40 and it should
    // be clamped to 1/CL=0.33 (and it is). However, we're getting a lot of
    // invalid combination printouts - that should NOT happen. It happens
    // because the algorithm generates in pp 0.441 which is > 1/CL, and that's
    // A BUG.
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.77},
        {1, 0.80},
        {2, 0.30},
        {3, 0.32},
    };
    int CL = 3;
#elif 0
    // BUG. We get the wrong percentages, and last node gets unusually
    // more. Is this related to the clamping of the incoming p's? probably
    // not.
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.66},
        {1, 0.66},
        {2, 0.30},
        {3, 0.30},
    };
    int CL = 3;
#elif 1
    // BUG. The error level is very low, but still we get a lot of invalid
    // combinations.
    // This happens because we get in pp a probability higher than 0.33 -
    // the first stage divided p node 1's causing me (which has the highest
    // deficit) to partitipate in the division twice and give it pp[1] = 0.12
    // and then later at the end we 2 remaining mixed nodes with surplus
    // 0.316, and add that to the 0.12 and get over 1/CL = 0.333....
    // Perhaps like we fixed the one-mixed-node-remaining case we also need
    // to fix this case? But it will be very messy to track for the different
    // me without running the full algorithm :-(
    std::vector<std::pair<int,float>> node_hit_rate {
        {0, 0.66},
        {1, 0.66},
        {2, 0.34},
        {3, 0.32},
    };
    int CL = 3;
#endif
// TODO: more tests - test also CL=1, CL=4, etc.
    int N = node_hit_rate.size();

    std::cout << "N=" << node_hit_rate.size() << " nodes, given hit rates:\n";
    for (auto& e : node_hit_rate) {
        std::cout << "    " << e.first << " " << e.second << "\n";
    }
    std::cout << "Asking for CL=" << CL << "\n\n";

    static thread_local std::default_random_engine random_engine;
    static thread_local std::uniform_int_distribution<> rand_node = std::uniform_int_distribution<>(0, N-1);

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
    // TODO: also count the times that a particular coordinator sent to
    // a particular CL combination 


    // debugging test - walk through one call for coordinator 2
//    miss_equalizing_combination(node_hit_rate, 2, CL);
//    exit(0);

    for (int i = 0; i < 1000000; i++) {
        // pick coordinator randomly, and then see which 2 nodes it will choose
        int coord = rand_node(random_engine);
//        std::cout << "coordinator " << coord << ": ";
        std::vector<int> c;
        if (reorder) {
            auto gen = miss_equalizing_combination(do_reorder(node_hit_rate, coord), 0, CL);
            c = gen.get();
        } else {
            auto gen = miss_equalizing_combination(node_hit_rate, coord, CL);
            c = gen.get();
        }
//        std::cout << "combination: " << c << "\n";
        // sort just for nicer printout (more useful with uniq)
        std::sort(c.begin(), c.end());
        bool invalid = false;
        int prev = -1;
        for(auto& s : c) {
//            std::cout << s << " ";
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
            std::cout << "invalid combination: " << c << "\n";
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
    }
    std::cout << "\n";
    if (count_invalid) {
        std::cout << "ERROR: Found invalid combinations: " << count_invalid << "\n";
    }
}
