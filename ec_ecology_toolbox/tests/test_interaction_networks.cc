//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
//#include "catch.hpp"

#include "../source/interaction_networks.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

PYBIND11_MODULE(lexicase, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("LexicaseFitness", &LexicaseFitness<emp::vector<double>>, "The lexicase function");
    //m.def("LexicaseFitness", [](emp::vector<emp::vector<double>> pop){return LexicaseFitness(pop);}, "The lexicase function");
}

/*
using org_t = emp::vector<int>;
using fit_map_t = emp::vector<double>;

TEST_CASE("VectorProduct", "[helpers]") {
    emp::vector<int> v({5,2,8});
    CHECK(VectorProduct(v) == 80);
}

TEST_CASE("FindHighest", "[helpers]") {
    emp::vector<org_t> v({{1,2,3}, {2, 1, 3}, {1,3,2}});
    auto highest = FindHighest(v, 0);
    CHECK(highest.size() == 1);
    emp::vector<int> res = highest[0];

    for (size_t i = 0; i < res.size(); i++) {
        CHECK( res[i] == v[1][i]);
    }

    highest = FindHighest(v, 1);
    CHECK(highest.size() == 1);
    res = highest[0];

    for (size_t i = 0; i < res.size(); i++) {
        CHECK( res[i] == v[2][i]);
    }

    highest = FindHighest(v, 2);
    CHECK(highest.size() == 2);

    emp::vector<emp::vector<double> > v2({{1.1,2.1,3.1}, {2.1, 1.1, 3.1}, {1.1,3.1,2.1}});
    auto highest2 = FindHighest(v2, 0);
    CHECK(highest2.size() == 1);
    emp::vector<double> res2 = highest2[0];

    for (size_t i = 0; i < res2.size(); i++) {
        CHECK( res2[i] == Approx(v2[1][i]));
    }

    highest2 = FindHighest(v2, 0, 1.0);
    CHECK(highest2.size() == 3);

    highest2 = FindHighest(v2, 1, 1.0);
    CHECK(highest2.size() == 2);


}

TEST_CASE("FilterImposssible" , "[helpers]") {
    emp::vector<emp::vector<double> > v({{1.5,2.1,3.1}, {1.5, 1.1, 1.1}, {1.1,1.1,1.1}, {1.5,0,0}, {0,0,0}, {1,1,3.1}});
    emp::vector<int> axes({0,1,2});
    emp::vector<emp::vector<double> > v2 = v;
    FilterImpossible(v2, axes);
    CHECK(v2.size() == 4);

    v2 = v;
    FilterImpossible(v2, axes, .4);
    CHECK(v2.size() == 5);

    axes = {1};
    v2 = v;
    FilterImpossible(v2, axes);
    CHECK(v2.size() == 1);

    axes = {2};
    v2 = v;
    FilterImpossible(v2, axes);
    CHECK(v2.size() == 2);
}

TEST_CASE("PruneAxes", "[helpers]") {
    emp::vector<org_t> v({{1,2,3}, {2, 1, 3}, {1,3,3}});
    emp::vector<int> axes({0,1,2});
    emp::vector<int> result(axes);
    PruneAxes(result, v);

    CHECK(result.size() == 2);
    CHECK(emp::Has(result, 0));
    CHECK(emp::Has(result, 1));

    emp::vector<emp::vector<double> > v2({{1.5,2.1,3.1}, {1.5, 1.1, 1.1}, {1.1,1.1,1.1}, {1.5,0,0}, {1.1,0,0}, {1.2,1,3.1}});
    result = axes;
    PruneAxes(result, v2);
    CHECK(result.size() == 3);
    result = axes;
    PruneAxes(result, v2, .4);
    CHECK(result.size() == 2);
    CHECK(emp::Has(result, 1));
    CHECK(emp::Has(result, 2));

}


TEST_CASE("Lexicase", "[selection_schemes]") {
    emp::vector<org_t> pop({{3,0,0}, {0, 3, 0}, {0,0,3}});
    fit_map_t fits = LexicaseFitness(pop);
    for (auto & o : fits) {
        CHECK(o == Approx(.33333333333));
    }

    pop = emp::vector<org_t>({{3,3,3}, {0, 1, 2}, {2,1,1}});
    fits = LexicaseFitness(pop);
    CHECK(fits[0] == 1);
    CHECK(fits[1] == 0);
    CHECK(fits[2] == 0);

    pop = emp::vector<org_t>({{3,3,3}, {3, 1, 2}, {2,1,1}});
    fits = LexicaseFitness(pop);
    CHECK(fits[0] == 1);
    CHECK(fits[1] == 0);
    CHECK(fits[2] == 0);

    pop = emp::vector<org_t>({{3,3,3}, {3, 3, 3}, {2,1,1}});
    fits = LexicaseFitness(pop);
    CHECK(fits[0] == Approx(.5));
    CHECK(fits[2] == 0);

    pop = emp::vector<org_t>({{3,1,2}, {1, 1, 2}, {2,1,1}});
    fits = LexicaseFitness(pop);
    CHECK(fits[0] == Approx(1));
    CHECK(fits[1] == 0);
    CHECK(fits[2] == 0);

    pop = emp::vector<org_t>({{3,1,2}, {1, 3, 2}, {2,3,1}});
    fits = LexicaseFitness(pop);
    CHECK(fits[0] == Approx(.5));
    CHECK(fits[1] == Approx(.333333));
    CHECK(fits[2] == Approx(.16666667));

    pop = emp::vector<org_t>({{3,1,2,1,1}, {1, 3, 2,1,1}, {2,3,1,1,1}});
    fits = LexicaseFitness(pop);
    CHECK(fits[0] == Approx(.5));
    CHECK(fits[1] == Approx(.333333));
    CHECK(fits[2] == Approx(.16666667));

    emp::vector<emp::vector<double>> pop_d = emp::vector<emp::vector<double>>({{3.1,1.1,2.1,1.1,1.1}, {1.1, 3.1, 2.1,1.1,1.1}, {2.1,3.1,1.1,1.1,1.1}});
    emp::vector<double> fits_d = LexicaseFitness(pop_d);
    CHECK(fits_d[0] == Approx(.5));
    CHECK(fits_d[1] == Approx(.333333));
    CHECK(fits_d[2] == Approx(.16666667));

    fits_d = LexicaseFitness(pop_d, emp::tools::Merge(Epsilon(1.0), DEFAULT));
    CHECK(fits_d[0] == Approx(0));
    CHECK(fits_d[1] == Approx(.25));
    CHECK(fits_d[2] == Approx(.75));


    emp::Random r;

    pop.resize(20);

    for (int org = 0; org < 20; ++org) {
        pop[org].resize(10);
        for (int loc = 0; loc < 10; ++loc) {
            pop[org].push_back(r.GetRandGeometric(.5));
        }
    } 

    fits = LexicaseFitness(pop);
    double total = 0;
    for (auto p : fits) {
        total += p;
    }
    CHECK(total == Approx(1));

    for (int org = 0; org < 20; ++org) {
        pop[org].resize(12);
        for (int loc = 0; loc < 12; ++loc) {
            pop[org].push_back(r.GetRandGeometric(.5));
        }
    } 

    fits = LexicaseFitness(pop);
    total = 0;
    for (auto p : fits) {
        total += p;
    }
    CHECK(total == Approx(1));

    for (int org = 0; org < 20; ++org) {
        pop[org].resize(13);
        for (int loc = 0; loc < 13; ++loc) {
            pop[org].push_back(r.GetRandGeometric(.5));
        }
    } 

    fits = LexicaseFitness(pop);
    total = 0;
    for (auto p : fits) {
        total += p;
    }
    CHECK(total == Approx(1));

    for (int org = 0; org < 20; ++org) {
        pop[org].resize(14);
        for (int loc = 0; loc < 14; ++loc) {
            pop[org].push_back(r.GetRandGeometric(.5));
        }
    } 

    fits = LexicaseFitness(pop);
    total = 0;
    for (auto p : fits) {
        total += p;
    }
    CHECK(total == Approx(1));

    for (int org = 0; org < 20; ++org) {
        pop[org].resize(18);
        for (int loc = 0; loc < 18; ++loc) {
            pop[org].push_back(r.GetRandGeometric(.5));
        }
    } 

    fits = LexicaseFitness(pop);
    total = 0;
    for (auto p : fits) {
        total += p;
    }
    CHECK(total == Approx(1));

}

TEST_CASE("Tournament", "[selection_schemes]") {
    emp::vector<org_t> pop = emp::vector<org_t>({{3,3,3}, {3, 3, 3}, {3,3,3}});
    emp::vector<double> fits = TournamentFitness(pop,  emp::tools::Merge(TournamentSize(2), DEFAULT));
    CHECK(fits[0] == Approx(.3333333));    

    fits = TournamentFitness(pop,  emp::tools::Merge(TournamentSize(3), DEFAULT));
    CHECK(fits[1] == Approx(.3333333));

    fits = TournamentFitness(pop,  emp::tools::Merge(TournamentSize(1), DEFAULT));
    CHECK(fits[2] == Approx(.3333333));


    pop = emp::vector<org_t>({{0,3,0}, {0, 3, 3}, {3,3,0}});
    fits = TournamentFitness(pop,  emp::tools::Merge(TournamentSize(3), DEFAULT));
    CHECK(fits[1] == Approx(13.0/27.0));    
    CHECK(fits[2] == Approx(13.0/27.0));    
    CHECK(fits[0] == Approx(1.0/27.0));    

    pop = emp::vector<org_t>({{0,3,0}, {0, 3, 3}, {3,3,3}});
    fits = TournamentFitness(pop,  emp::tools::Merge(TournamentSize(3), DEFAULT));
    CHECK(fits[1] == Approx(7.0/27.0));    
    CHECK(fits[2] == Approx(19.0/27.0));    
    CHECK(fits[0] == Approx(1.0/27.0));    

    pop = emp::vector<org_t>({{0,3,3}, {3, 3, 3}, {3,3,3}});
    fits = TournamentFitness(pop,  emp::tools::Merge(TournamentSize(3), DEFAULT));
    CHECK(fits[1] == Approx(13.0/27.0));    
    CHECK(fits[0] == Approx(1.0/27.0));    

    fits = TournamentFitness(pop,  emp::tools::Merge(TournamentSize(2), DEFAULT));
    CHECK(fits[2] == Approx(.44444));    
    CHECK(fits[0] == Approx(.11111));    

    fits = TournamentFitness(pop,  emp::tools::Merge(TournamentSize(1), DEFAULT));
    CHECK(fits[2] == Approx(.3333333));    
    CHECK(fits[0] == Approx(.3333333));    

}

TEST_CASE("Fitness sharing", "[selection_schemes]") {
    emp::vector<org_t> pop = emp::vector<org_t>({{3,3,3}, {3, 3, 3}, {3,3,3}});
    all_attrs settings = DEFAULT;
    fit_map_t fits = SharingFitness(pop, settings);
    CHECK(fits[0] == Approx(.3333333));

    pop = emp::vector<org_t>({{3,1,2,1,1}, {1, 3, 2,1,1}, {2,3,1,1,1}});
    fits = SharingFitness(pop, DEFAULT);
    CHECK(fits[0] == Approx(.5555556));
    CHECK(fits[1] == Approx(.333333));
    CHECK(fits[2] == Approx(.111111));

    pop = emp::vector<org_t>({{10,1}, {1, 10}, {1,1}});
    fits = SharingFitness(pop, DEFAULT);
    CHECK(fits[0] == Approx(.444444));
    CHECK(fits[1] == Approx(.444444));
    CHECK(fits[2] == Approx(.111111));

    emp::Random r;

    pop.resize(20);
    for (int org = 0; org < 20; ++org) {
        pop[org].resize(10);
        for (int loc = 0; loc < 10; ++loc) {
            pop[org].push_back(r.GetRandGeometric(.5));
        }
    } 

    fits = SharingFitness(pop, DEFAULT);
    double total = 0;
    for (auto p : fits) {
        total += p;
    }
    CHECK(total == Approx(1));

}

TEST_CASE("Eco-EA", "[selection_schemes]") {
    emp::vector<org_t> pop = emp::vector<org_t>({{3,3,3}, {3, 3, 3}, {3,3,3}});
    all_attrs settings = DEFAULT;
    fit_map_t fits = EcoEAFitness(pop, settings);
    CHECK(fits[1] == Approx(.3333333));

    pop = emp::vector<org_t>({{3,1,2,1,1}, {1, 3, 2,1,1}, {2,3,1,1,1}});
    fits = EcoEAFitness(pop, DEFAULT);
    CHECK(fits[0] == Approx(5.0/9.0));
    CHECK(fits[1] == Approx(2.0/9.0));
    CHECK(fits[2] == Approx(2.0/9.0));

    pop = emp::vector<org_t>({{10,1,2,1,1}, {1, 3, 2,1,1}, {2,3,1,1,1}, {2,1,1,1,1}});
    fits = EcoEAFitness(pop, DEFAULT);
    CHECK(fits[0] == Approx(7.0/16.0));
    CHECK(fits[1] == Approx(2.0/16.0));
    CHECK(fits[2] == Approx(2.0/16.0));
    CHECK(fits[3] == Approx(5.0/16.0));

    emp::Random r;
    pop.resize(20);
    for (int org = 0; org < 20; ++org) {
        pop[org].resize(10);
        for (int loc = 0; loc < 10; ++loc) {
            pop[org].push_back(r.GetRandGeometric(.5));
        }
    } 

    fits = EcoEAFitness(pop, DEFAULT);
    double total = 0;
    for (auto p : fits) {
        total += p;
    }
    CHECK(total == Approx(1));

}

TEST_CASE("Calc competition", "[helpers]") {
    emp::vector<org_t> pop = emp::vector<org_t>({{1,3,1}, {3, 1, 1}, {1,1,3}});
    // all_attrs settings = DEFAULT;
    std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> test_fun = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {
        fit_map_t base_fit_map;
        for (int i = 0; i < pop.size(); i++) {
            base_fit_map.push_back(1.0);
        }
        return base_fit_map;
    };

    emp::WeightedGraph g = CalcCompetition(pop, test_fun);
    auto weights = g.GetWeights();
    for (auto vec : weights) {
        for (auto val : vec) {
            CHECK(val == 0); 
        }
    }

    test_fun = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {
        fit_map_t base_fit_map;
        for (int i = 0; i < pop.size(); i++) {
            base_fit_map.push_back(1.0);
        }

        if (emp::Has(pop, {1,3,1})) {
            base_fit_map[1] = 0;
        }
        return base_fit_map;
    };

    g = CalcCompetition(pop, test_fun);
    CHECK(g.GetWeight(0,1) == -1);

    pop = emp::vector<org_t>({{1,0,1}, {0,1,0}, {0,0,1}, {2,0,0}});

    g = CalcCompetition(pop, do_lexicase<org_t>);
    CHECK(g.GetWeight(0,3) == 0);
    CHECK(g.GetWeight(3,0) == Approx(-1.0/3.0));
    CHECK(g.GetWeight(0,2) == Approx(-1.0/3.0));
    CHECK(g.GetWeight(2,0) == 0);

    g = CalcCompetition(pop, do_tournament<org_t>);

    CHECK(g.GetWeight(1,3) == 0);

    // The following are regression tests (I have not done the math to confirm the numbers are right
    // but they seem plausible and accuracy is being tested elsewhere)

    g = CalcCompetition(pop, do_eco_ea<org_t>, emp::tools::Merge(NicheWidth(1.0), DEFAULT));

    CHECK(g.GetWeight(0,0) == 0);
    CHECK(g.GetWeight(0,1) == Approx(.1875));
    CHECK(g.GetWeight(0,2) == Approx(.0625));
    CHECK(g.GetWeight(0,3) == Approx(.125));
    CHECK(g.GetWeight(1,0) == 0);
    CHECK(g.GetWeight(1,1) == 0);
    CHECK(g.GetWeight(1,2) == 0);
    CHECK(g.GetWeight(1,3) == Approx(.125));
    CHECK(g.GetWeight(2,0) == 0);
    CHECK(g.GetWeight(2,1) == Approx(.125));
    CHECK(g.GetWeight(2,2) == 0);
    CHECK(g.GetWeight(2,3) == Approx(.125));
    CHECK(g.GetWeight(3,0) == 0);
    CHECK(g.GetWeight(3,1) == 0);
    CHECK(g.GetWeight(3,2) == 0);
    CHECK(g.GetWeight(3,3) == 0);

    g = CalcCompetition(pop, do_sharing<org_t>);

    CHECK(g.GetWeight(0,0) == 0);
    CHECK(g.GetWeight(0,1) == Approx(-.0625));
    CHECK(g.GetWeight(0,2) == Approx(-.1875));
    CHECK(g.GetWeight(0,3) == 0);
    CHECK(g.GetWeight(1,0) == 0);
    CHECK(g.GetWeight(1,1) == 0);
    CHECK(g.GetWeight(1,2) == Approx(-.125));
    CHECK(g.GetWeight(1,3) == 0);
    CHECK(g.GetWeight(2,0) == 0);
    CHECK(g.GetWeight(2,1) == 0);
    CHECK(g.GetWeight(2,2) == 0);
    CHECK(g.GetWeight(2,3) == 0);
    CHECK(g.GetWeight(3,0) == Approx(-.125));
    CHECK(g.GetWeight(3,1) == Approx(-.125));
    CHECK(g.GetWeight(3,2) == Approx(-.125));
    CHECK(g.GetWeight(3,3) == 0);

}
*/