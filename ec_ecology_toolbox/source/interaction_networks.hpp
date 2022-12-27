#ifndef _INTERACTION_NETWORKS_H_
#define _INTERACTION_NETWORKS_H_

//#include <pybind11/pybind11.h>

#include <functional>
#include <map>
#include <limits>
#include <unordered_set>

#include "base/vector.hpp"
#include "base/assert.hpp"
#include "math/Random.hpp"
#include "math/random_utils.hpp"
#include "datastructs/vector_utils.hpp"
#include "datastructs/map_utils.hpp"
#include "tools/string_utils.hpp"
#include "datastructs/set_utils.hpp"
#include "math/math.hpp"
#include "math/distances.hpp"
#include "tools/attrs.hpp"
#include "datastructs/Graph.hpp"
#include "Evolve/Resource.hpp"

#include <math.h>  



DEFINE_ATTR(SigmaShare);
DEFINE_ATTR(Alpha);
DEFINE_ATTR(Cost);
DEFINE_ATTR(Cf);
DEFINE_ATTR(NicheWidth);
DEFINE_ATTR(MaxScore);
DEFINE_ATTR(ResourceInflow);
DEFINE_ATTR(ResourceOutflow);
DEFINE_ATTR(MaxBonus);
DEFINE_ATTR(TournamentSize);
DEFINE_ATTR(Epsilon);

constexpr auto DEFAULT{MakeAttrs(SigmaShare(8.0),
                                 Alpha(1.0),
                                 Cost(1.0),
                                 Cf(.0025),                                 
                                 NicheWidth(3.0),
                                 MaxScore(10.0),
                                 ResourceInflow(2000.0),
                                 ResourceOutflow(.01),
                                 MaxBonus(5.0),
                                 TournamentSize(2),
                                 Epsilon(0.0))};

using all_attrs = emp::tools::Attrs<typename SigmaShare::value_t<double>, typename Alpha::value_t<double>, 
                        typename Cost::value_t<double>, typename Cf::value_t<double>,
                        typename NicheWidth::value_t<double>, typename MaxScore::value_t<double>,
                        typename ResourceInflow::value_t<double>, typename ResourceOutflow::value_t<double>,
                        typename MaxBonus::value_t<double>, typename TournamentSize::value_t<int>, typename Epsilon::value_t<double> >;


// from https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T, class U>
typename std::enable_if<!std::numeric_limits<T>::is_integer || !std::numeric_limits<U>::is_integer, bool>::type
    almost_equal(T x, U y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}

// integer version - just use the equals operator
template<class T, class U>
typename std::enable_if<std::numeric_limits<T>::is_integer && std::numeric_limits<U>::is_integer, bool>::type
    almost_equal(T x, U y, int ulp)
{
    return x == y;
}

/// Find the elements in @param pop with the highest value on
/// @param axis. PHEN_T should be a vector of some numeric type.
/// If epsilon lexicase is being used, specify a value of @param epsilon
/// to include everything within that range of the highest value
template <typename PHEN_T>
emp::vector<int> FindHighestIndices(emp::vector<PHEN_T> & pop, int axis, double epsilon = 0) {
    double best = std::numeric_limits<double>::lowest();
    emp::vector<int> winners;

    for (size_t i = 0; i < pop.size(); i++) {
        if (pop[i][axis] > best) {
            best = pop[i][axis];
            winners.resize(0);
            winners.push_back(i);
        } else if (almost_equal(pop[i][axis], best, 5)) {
            winners.push_back(i);
        }
    }

    if (epsilon) {
        winners.resize(0);
        for (size_t i = 0; i < pop.size(); i++) {
            if (pop[i][axis] + epsilon >= best) {
                winners.push_back(i);
            }
        }
    }
    return winners;
}

template <typename PHEN_T>
emp::vector<PHEN_T> FindHighest(emp::vector<PHEN_T> & pop, int axis, double epsilon = 0) {
    // From https://stackoverflow.com/questions/9650679/c-select-a-subset-of-a-stdvector-based-predefined-element-indices

    emp::vector<int> highest = FindHighestIndices(pop, axis, epsilon);
    std::vector<PHEN_T> result(highest.size());

    std::transform(highest.begin(), highest.end(), result.begin(), [&pop](size_t pos) {return pop[pos];});
    return result;
}

/// Multiply all of the ints in @param nums together.
long int VectorProduct(const emp::vector<int> & nums) {
    long int result = 1;

    for (int num : nums) {
        result *= num;
    }

    return result;
}

template <typename PHEN_T>
void FilterImpossible(emp::vector<PHEN_T> & pop, emp::vector<int> & axes, double epsilon = 0) { 
    emp::vector<int> best;
    emp::vector<int> temp;
    for (int ax : axes) {
        temp = FindHighestIndices(pop, ax, epsilon);
        best.reserve(best.size() + distance(temp.begin(),temp.end()));
        best.insert(best.end(),temp.begin(),temp.end());       
    }

    best = emp::RemoveDuplicates(best);
    std::vector<PHEN_T> result(best.size());

    std::transform(best.begin(), best.end(), result.begin(), [&pop](size_t pos) {return pop[pos];});

    pop = result;
}

template <typename PHEN_T>
void PruneAxes(emp::vector<int> & axes, emp::vector<PHEN_T> & pop, double epsilon = 0) {
    
    emp::vector<int> to_remove;

    for (int ax : axes) {
        // std::cout << "Evaluating " << ax << std::endl;
        bool include = false;
        double best = pop[0][ax];
        double lowest = pop[0][ax];
        for (PHEN_T & org : pop) {
            if (org[ax] > best) {
                // std::cout << org[ax] << " > "<< best << std::endl;
                if (org[ax] > lowest + epsilon) {
                    // std::cout <<" Including  "<< ax << std::endl;
                    include = true;
                    break;
                }
                best = org[ax];
            } else if (org[ax] < lowest) {
                // std::cout << org[ax] << " < "<< lowest << std::endl;
                if (org[ax] + epsilon < best) {
                    // std::cout <<" Including  "<< ax << std::endl;
                    include = true;
                    break;
                }
                lowest = org[ax];
            }
        }

        if (!include) {
            to_remove.push_back(ax);
        }
    }
    for (int ax : to_remove) {
        axes.erase(std::remove(axes.begin(), axes.end(), ax), axes.end());
    }
}

template <typename PHEN_T>
void HandleTwoOrgs(std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> & winners, emp::vector<int> axes, emp::vector<int> perm_levels, double epsilon = 0, emp::vector<int> active_set = {}) {
    double wins = 0;
    double losses = 0;
    for (int ax : axes) {
        if (winners[0][ax] > winners[1][ax] + epsilon) {
            wins++;
        } else if (winners[1][ax] > winners[0][ax] + epsilon) {
            losses++;
        }
    }
    if (wins > 0 || losses > 0) {
        fit_map[winners[0]] += (wins/(wins+losses)) * (1.0/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
        fit_map[winners[1]] += (losses/(wins+losses)) * (1.0/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
    } else {
        fit_map[winners[0]] += .5 * (1.0/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
        fit_map[winners[1]] += .5 * (1.0/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);                
    }
}

// template <typename PHEN_T>
// double SizeTieSet(PHEN_T & org, std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> winners, emp::vector<int> axes, emp::vector<int> perm_levels, double epsilon = 0) {
//     emp::vector<int> temp = FindHighestIndices(winners, axes[0], epsilon);
//     std::set<int> intersect_set(temp.begin(), temp.end());
//     for (int ax : axes) {
//         intersect_set = emp::intersection(intersect_set, FindHighestIndices(winners, ax, epsilon));
//     }

//     return (double)intersect_set.size();
// }

// template <typename PHEN_T>
// void OrgWinCondition(PHEN_T & org, std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> winners, emp::vector<int> axes, emp::vector<int> perm_levels, double epsilon = 0) {
//     // if (axes.size() == 1) {
//     //     fit_map[org] += (1.0/VectorProduct(perm_levels)) * (1.0/(double)winners.size());
//     //     return;
//     // }
    
//     // PruneAxes(axes, winners, epsilon);

//     emp::vector<int> wins;
//     emp::vector<int> ties;
//     emp::vector<int> losses;
//     for (int ax : axes) {
//         bool win = true;
//         bool tie = true;
//         for (auto & other : winners) {
//             if (other[ax] > org[ax]+epsilon) {
//                 win = false;
//                 tie = false;
//                 break;
//             } else if (!(org[ax] > other[ax] + epsilon) && (other != org)) {
//                 win = false;
//             }
//         }
//         if (win) {
//             wins.push_back(ax);
//         } else if(tie) {
//             ties.push_back(ax);
//         } else {
//             losses.push_back(ax);
//         }
//     }


//     if (wins == axes.size()) {
//         fit_map[org] += 1.0/VectorProduct(perm_levels);
//     } else {
//         perm_levels.push_back(axes.size());
//         fit_map[org] += (wins/(double)axes.size()) * 1.0/VectorProduct(perm_levels);
//         // double tie_prob = 0.0;
//         // for (size_t i = 0; i < ties.size()-1; i++) {
//         //     tie_prob += (wins/(double)axes.size()) * pow(losses/(double)axes.size(), i);
//         for (int ax : ties) {
//             std::cout<< "exploring " << ax << emp::to_string(ties) << std::endl; 
//             emp::vector<int> next_axes = axes;
//             next_axes.erase(std::remove(next_axes.begin(), next_axes.end(), ax), next_axes.end());
//             OrgWinCondition(org, fit_map, FindHighest(winners, ax, epsilon), next_axes, perm_levels, epsilon);
//         }
//         // }
//         // tie_prob += pow(losses/(double)axes.size(), ties.size()-1) * 1.0/SizeTieSet(org, fit_map, winners, ties, perm_levels, epsilon); // Actual tie
//         // fit_map[org] += tie_prob * 1.0/VectorProduct(perm_levels);
//     }

// }

// template <typename PHEN_T>
// void HandleThreeOrgs(std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> & winners, emp::vector<int> axes, emp::vector<int> perm_levels, double epsilon = 0, emp::vector<int> active_set = {}) {
//     for (PHEN_T & org : winners) {
//         OrgWinCondition(org, fit_map, winners, axes, perm_levels, epsilon);
//     }

// }

template <typename PHEN_T>
void TraverseDecisionTree(std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> & pop, emp::vector<int> axes, emp::vector<int> perm_levels, double epsilon = 0, emp::vector<int> active_set = {}) {
    // std::cout << "begining round of recursion " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    // std::cout << emp::to_string(pop) << emp::to_string(axes) << std::endl;
 
    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org]+=1.0/((double)winners.size()*VectorProduct(perm_levels));
        }
        return;
    }

    if (pop.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org]+=1.0/VectorProduct(perm_levels);
        }
        return;
    }

    PruneAxes(axes, pop, epsilon);
    perm_levels.push_back(axes.size());

    FilterImpossible(pop, axes, epsilon);

    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org]+=1.0/((double)winners.size()*VectorProduct(perm_levels));
        }
        return;
    }

    // std::cout << "post processing: " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    for (int ax : axes) {
        // if (perm_levels.size() == 1) {
        //     std::cout << "Axis: " << ax << " out of " << emp::to_string(axes) << std::endl;
        // }
        emp::vector<PHEN_T> winners = FindHighest(pop, ax, epsilon);
        emp::vector<int> next_axes = axes;
        next_axes.erase(std::remove(next_axes.begin(), next_axes.end(), ax), next_axes.end());
        FilterImpossible(winners, next_axes, epsilon);

        if (winners.size() == 1) { // Not a tie
            // std::cout << "1 winner: " << emp::to_string(winners[0]) << " Controls " << (double)emp::Factorial(axes.size() - 1)<< std::endl;
            fit_map[winners[0]] += 1.0/VectorProduct(perm_levels);//(double)emp::Factorial(axes.size() - 1);
        } else if (winners.size() == 2) { // optimization
            HandleTwoOrgs(fit_map, winners, next_axes, perm_levels, epsilon);
        // } else if (winners.size() < 100 && next_axes.size() > 2) { // optimization
        //     HandleThreeOrgs(fit_map, winners, next_axes, perm_levels, epsilon);
        } else { // tie
            TraverseDecisionTree(fit_map, winners, next_axes, perm_levels, epsilon);
        }
    }
}

template <typename PHEN_T>
emp::vector<double> LexicaseFitness(emp::vector<PHEN_T> & pop) {

    emp_assert(pop.size() > 0);
    std::map<PHEN_T, double> fit_map;
    size_t n_funs = pop[0].size();

    for (PHEN_T & org : pop) {
        fit_map[org] = 0.0;
    }

    emp::vector<PHEN_T> de_dup_pop = emp::RemoveDuplicates(pop);
    TraverseDecisionTree(fit_map, de_dup_pop, emp::NRange(0, (int)n_funs), {});

    for (PHEN_T & org : de_dup_pop) {
        fit_map[org] /= emp::Count(pop, org);
        // std::cout << "Pre div: " << fit_map[org] << std::endl;
        // fit_map[org] /= emp::Factorial(n_funs); // convert to proportion of "islands"
        // std::cout << "Post div: " << fit_map[org] << " Divided by: " << (double)emp::Factorial(n_funs) << std::endl;
    }

    emp::vector<double> result;
    for (PHEN_T & org : pop) {
        result.push_back(fit_map[org]);
    }

    return result;
}

void TournamentHelper(emp::vector<double> & fit_map, int t_size = 2){

    emp::vector<double> base_fit_map = fit_map;
    double pop_size = base_fit_map.size();

    for (size_t i = 0; i < fit_map.size(); i++) {
        double less = 0.0;
        double equal = 0.0; 

        for (auto & org2 : base_fit_map) {
            if (almost_equal(org2, base_fit_map[i], 10)) {
                equal++;
            } else if (org2 < base_fit_map[i]) {
                less++;
            }
        }

        long double p_less = less/pop_size;
        long double p_equal = equal/pop_size;
        long double p_self = 1/pop_size;

        fit_map[i] = (pow(p_equal + p_less, t_size) - pow(p_less, t_size)) * p_self/p_equal;
    }
}

template <typename PHEN_T>
emp::vector<double> RandomFitness(emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT) {
    emp_assert(pop.size() > 0);
    emp::vector<double> fit_map;
    for (PHEN_T & org : pop) {
        fit_map.push_back(1.0/(double)pop.size());
    }
    TournamentHelper(fit_map, TournamentSize::Get(attrs));
    return fit_map;
}

template <typename PHEN_T>
emp::vector<double> TournamentFitness(emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT) {
    emp_assert(pop.size() > 0);
    emp::vector<double> fit_map;
    for (PHEN_T & org : pop) {
        fit_map.push_back(emp::Sum(org));
    }
    TournamentHelper(fit_map, TournamentSize::Get(attrs));
    return fit_map;
}

template <typename PHEN_T>
emp::vector<double> SharingFitness(emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT) {
    // std::cout << "SHARING" << std::endl;
    emp::vector<double> fit_map;

    // std::cout << SigmaShare::Get(attrs) << std::endl;

    for (size_t i = 0; i < pop.size(); i++) {
        fit_map.push_back(1.0);
    }

    for (size_t i = 0; i < pop.size(); i++) {
        double niche_count = 0;

        for (PHEN_T & org2 : pop) {
            // Sharing function is euclidean distance
            // we could make this configurable
            double dist = emp::EuclideanDistance(pop[i], org2);
            if (dist < SigmaShare::Get(attrs)) {
                niche_count += 1 - pow((dist/SigmaShare::Get(attrs)), Alpha::Get(attrs));
            } 
        }

        // Don't worry, niche_count will always be at least one because each individual
        // increases its own niche count by 1
        fit_map[i] = emp::Sum(pop[i]) / niche_count;
    }
    TournamentHelper(fit_map, TournamentSize::Get(attrs));

    return fit_map;    
};

template <typename PHEN_T>
emp::vector<double> EcoEAFitness(emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT) {
    // std::cout << "ECO-EA" << std::endl;
    emp_assert(pop.size() > 0);

    emp::vector<double> fit_map;

    size_t n_funs = pop[0].size();

    for (size_t i = 0; i < pop.size(); i++) {
        fit_map.push_back(1.0);
    }

    for (int axis : emp::NRange(0, (int)n_funs)) {
        double res = ResourceInflow::Get(attrs)/(double)pop.size(); // Average resource available per individual
        double count = 0;

        for (PHEN_T & org : pop) {
            if (org[axis] >= NicheWidth::Get(attrs)) {
                count+= std::min(Cf::Get(attrs)*res*pow(org[axis]/MaxScore::Get(attrs),2.0) - Cost::Get(attrs), MaxBonus::Get(attrs));
            }
        }

        count /= (double) pop.size(); // average resource use per pop member
        res -= count; // how much of average available resources are used on average
        res = std::max(res, 0.0); // Can't have a negative amount of a resource

        for (size_t i = 0; i < pop.size(); i++) {
            if (pop[i][axis] >= NicheWidth::Get(attrs)) {
                fit_map[i] *= pow(2,std::min(Cf::Get(attrs)*res*pow(pop[i][axis]/MaxScore::Get(attrs),2.0) - Cost::Get(attrs), MaxBonus::Get(attrs)));
            }   
        }
    }

    TournamentHelper(fit_map, TournamentSize::Get(attrs));

    return fit_map;
};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_lexicase = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return LexicaseFitness(pop,attrs);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_tournament = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return TournamentFitness(pop,attrs);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_eco_ea = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return EcoEAFitness(pop,attrs);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_sharing = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return SharingFitness(pop,attrs);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_random = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return RandomFitness(pop,attrs);};


template <typename PHEN_T>
emp::WeightedGraph CalcCompetition(emp::vector<PHEN_T> pop, 
                        std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> fit_fun,
                        all_attrs attrs=DEFAULT) {

    emp::WeightedGraph effects(pop.size());


    emp::vector<double> fitnesses = fit_fun(pop, attrs);

    for (size_t i = 0; i < pop.size(); i++) {
        effects.SetLabel(i, emp::to_string(pop[i]));
        // std::cout << effects.GetLabel(i) << std::endl;

        emp::vector<PHEN_T> curr = pop;
        for (size_t ax = 0; ax < curr[i].size(); ax++) {
            curr[i][ax] = 0; // Replace org with null org so pop size stays same
        }

        emp::vector<double> new_fits = fit_fun(curr, attrs);
        for (size_t j = 0; j < pop.size(); j++ ) {
            if (i == j) {continue;}

            // In terms of floating-point imprecision issues, it's much better
            // to check for equality than doing the subtraction and checking
            // for the result to be 0
            if (!almost_equal(fitnesses[j], new_fits[j], 10)) {
                effects.AddEdge(i, j, fitnesses[j] - new_fits[j]);
            }
            // std::cout << effect << std::endl;
        }
    }

    return effects;
}



#endif
