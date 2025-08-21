// g++ -o checker Checker.cpp -fopenmp -lgmp -lgmpxx -O3

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <functional>
#include <numeric>
#include <utility>
#include <cstdlib>
#include <gmpxx.h>
#include <omp.h>
#include "progressbar.hpp"

typedef std::vector<mpz_class> Residues;
typedef std::vector<std::pair<int, std::set<int>>> Filters;

namespace IO
{
    Residues readResidues(const std::string& filepath)
    {
        std::ifstream fin(filepath);
        Residues R;
        mpz_class r;

        while(fin >> r)
            R.push_back(r);

        fin.close();
        return R;
    }

    Filters readFilters(const std::string& filepath,
                        const int& filterCount = 140000)
    {
        std::ifstream fin(filepath);
        Filters F;
        int p = -1, n = 0, filtersRead = 0;

        while(fin >> n)
        {
            if(p == -1)
                F.push_back(std::pair<int, std::set<int>>(p = n, std::set<int>()));
            else if(n == -1)
            {
                p = -1;
                if(++filtersRead >= filterCount)
                    break;
            }
            else
                F.back().second.insert(n);
        }

        fin.close();
        return F;
    }
}

// Edit this to wherever you have your residues and filters file.
const Residues R = IO::readResidues("resources/Residues.txt");
const Filters F = IO::readFilters("resources/Filters.txt");

// Edit this to your G if you're using other residues.
const mpz_class G = 25878772920;

namespace Conjecture
{
    [[maybe_unused]]
    bool validate(const mpz_class& p,
                  const mpz_class& x,
                  const mpz_class& y,
                  const mpz_class& z)
    {
        mpz_class lhs = 4 * x * y * z;
        mpz_class rhs = p * (x * y + y * z + z * x);
        return lhs == rhs;
    }

    bool filterRemove(const mpz_class& n)
    {
        for(const auto& [p, f] : F)
        {
            mpz_class r_large = n % p;
            int r = r_large.get_ui();
            if(f.count(r) != 0)
                return true;
        }

        return false;
    }

    bool checkIndividual(const mpz_class& n)
    {
        if(filterRemove(n))
            return true;

        // More checks can be added here if needed.

        return false;
    }

    std::vector<mpz_class> checkBatch(const mpz_class& k)
    {
        std::vector<mpz_class> hardResidues;

        for(const auto& r : R)
        {    
            mpz_class true_r = r + G * k;
            if(!checkIndividual(true_r))
                hardResidues.push_back(true_r);
        }

        return hardResidues;
    }
}

int main(int argc, char* argv[])
{
    /**
     * @note These are the batches that must be checked for 10^18.
     * int firstBatch = 3864170, lastBatch = 38641709; 
     **/

    std::ofstream fout(argv[1]);
    int firstBatch = strtol(argv[2], nullptr, 10), lastBatch = strtol(argv[3], nullptr, 10);
    progressbar progressBar(lastBatch - firstBatch);

    std::cout << "Checking the conjecture for batches " << firstBatch << '-' << lastBatch << '\n';

    #pragma omp parallel for
    for(int batch = firstBatch; batch <= lastBatch; batch++)
    {
        std::vector<mpz_class> misses = Conjecture::checkBatch(batch);

        #pragma omp critical
        {
            for(const auto& miss : misses)
                fout << miss << '\n';
                
            progressBar.update();
        }
    }

    return 0;
}
