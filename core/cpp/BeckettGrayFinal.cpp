using namespace std;

#include <fstream>  // File stream: used for reading from or writing to files on disk
#include <iostream> // Input/Output stream: standard console printing (cout, cerr)
#include <vector>   // Dynamic array container: provides resizable storage for data
#include <cmath>    // Mathematical functions: provides pow, sqrt, log, etc.
#include <cstdio>   // C-style Input/Output: used for fast I/O like printf or scanf
#include <cstdlib>  // Standard library: general purpose functions (memory, rand, exit)
#include <ctime>    // C-style Time: used for basic timing and seeding random numbers
#include <bitset>   // Bit manipulation: useful for fixed-size binary representations

#include "EvalPerf.hpp"

// --- Core Algorithm Variables ---
int n;       // Number of bits (dimension of the hypercube)
int n2;      // Total number of nodes (2^n)
int s;       // Current sequence length during recursion
int smax;    // Maximum sequence length reached so far (tracking progress)
int nbAppel; // Recursive function call counter (performance monitoring)

// --- Sequence Storage ---
vector<int> bgc;                  // Current Beckett Gray Code sequence being constructed
vector<vector<int>> listeCodeBGC; // Collection of all valid Beckett Gray Codes found

// --- State Tracking & Constraints ---
vector<bool> avail; // Visitation map (True if the binary integer is unvisited)
vector<int> old;    // FIFO-like buffer of bit indices currently set to 1
vector<int> oldS;   // Timestamps (step indices) when bits were flipped to 1
vector<int> nbMax;  // Balance constraints: max allowed consecutive 1s per bit position

// --- Precomputed Lookup Tables (Optimization) ---
vector<vector<int>> bi;    // Precomputed binary representation for each integer < 2^n
vector<vector<int>> flip;  // Adjacency matrix: flip[i][j] is the neighbor of i by toggling bit j
vector<int> nombreVoisins; // Current degree of each vertex (number of unvisited neighbors)
vector<int> V;             // Hamming weight: V[i] stores the number of set bits (1s) in i
vector<int> nbEltEns;      // Population count: number of remaining elements for each Hamming weight

// --- Performance & Statistics ---
EvalPerf PE; // Performance evaluation tool (timers and hardware counters)
int bmin;    // Minimum "max-balance" found across all generated codes
int bmax;    // Maximum "max-balance" found across all generated codes

/**
 * Converts an integer into a binary representation as a vector of bits.
 * @param x The integer to convert.
 * @param n The number of bits to extract (starting from the right).
 * @return A vector of n integers (0 or 1).
 */
vector<int> Binary(int x, int n)
{
    // Convert x to a 32-bit binary string representation
    string b = bitset<32>(static_cast<unsigned long long>(x)).to_string();

    vector<int> res;
    // Extract the last n bits by iterating from the (32-n)th character
    for (int i = b.size() - n; i < b.size(); i++)
    {
        // Convert char ('0' or '1') to int by subtracting the ASCII value of '0'
        res.push_back(b[i] - '0');
    }
    return res;
}

/**
 * Populates the Hamming weight vector (V) and counts elements per weight.
 * @param n The number of bits (dimension).
 * @return A vector where res[i] is the number of set bits (1s) in integer i.
 */
vector<int> listeDes1(int n)
{
    vector<int> res;

    for (int i = 0; i < pow(2, n); i++)
    {
        int cmpt = 0;
        int nb = i;

        // Brian Kernighan's-like bit counting algorithm
        while (nb > 0)
        {
            cmpt += nb & 1; // Check if the least significant bit is 1
            nb >>= 1;       // Right shift to process the next bit
        }
        res.push_back(cmpt);

        // Increment the population count for this specific Hamming weight
        nbEltEns[cmpt]++;
    }
    return res;
}

/**
 * Generates a lookup table containing binary representations of all integers up to 2^n - 1.
 * @param n The number of bits for each representation.
 * @return A 2D vector where res[i] is the binary vector (0s and 1s) of integer i.
 */
vector<vector<int>> listeBinaire(int n)
{
    vector<vector<int>> res;
    for (int i = 0; i < pow(2, n); i++)
    {
        // Compute and store the fixed-length bit vector for each integer
        res.push_back(Binary(i, n));
    }
    return res;
}

/**
 * Checks if a given integer is a power of two using bitwise logic.
 * @param n The integer to check.
 * @return True if n is a power of two, False otherwise.
 */
bool pw2(int n)
{
    // A power of two has exactly one bit set to 1.
    // The operation (n & (n - 1)) clears the lowest set bit.
    // If the result is 0 and n was not 0, then n is a power of two.
    return n != 0 && (n & (n - 1)) == 0;
}

/**
 * Generates an adjacency list for all integers in the 2^n space.
 * * For each integer i, it computes all neighbors reachable by flipping exactly one bit.
 * This effectively maps the edges of an n-dimensional hypercube.
 *
 * @param n The number of bits (dimension of the hypercube).
 * @return A 2D vector where res[i][j] is the result of flipping the j-th bit of i.
 */
vector<vector<int>> listeFlip(int n)
{
    vector<vector<int>> res;

    for (int i = 0; i < pow(2, n); i++)
    {
        res.push_back(vector<int>());
        for (int j = 0; j < n; j++)
        {
            // Bitwise XOR with a single-bit mask (1 << j) toggles the j-th bit
            res[i].push_back(i ^ (1 << j));
        }
    }
    return res;
}

/**
 * Prints a boolean vector to the standard output.
 * * @param l The vector of booleans to be displayed.
 */
void toStringBool(vector<bool> l)
{
    // Iterate through the vector and print 1 for True, 0 for False
    for (bool b : l)
    {
        cout << b << " ";
    }
    cout << endl; // New line after printing all elements
}

/**
 * Prints an integer vector to the standard output.
 * * @param l The vector of integers to be displayed.
 */
void toStringInt(vector<int> l)
{
    // Iterate through the vector and print each integer separated by a space
    for (int i : l)
    {
        cout << i << " ";
    }
    // Print a newline character at the end of the sequence
    cout << endl;
}

/**
 * Analyzes a specific Gray Code to find its maximum run-length of consecutive 1s.
 * Updates the global bmin and bmax if a new best or worst balance is found.
 * * @param i The index of the Gray Code to analyze within listeCodeBGC.
 */
void compteLongeur1(int i)
{
    int max = 0;

    // Check each bit position (column) k
    for (int k = 0; k < n; k++)
    {
        int cmpt = 0;

        // Traverse the entire sequence (n2 = 2^n)
        for (int j = 0; j < n2; j++)
        {
            // Access the k-th bit of the j-th number in the sequence
            if (bi[listeCodeBGC[i][j]][k] == 1)
            {
                cmpt++;
            }
            else
            {
                // Sequence of 1s broken, update the max if necessary
                if (cmpt > max)
                {
                    max = cmpt;
                }
                cmpt = 0;
            }
        }
    }

    // Update global maximum balance if this code has a longer run of 1s
    if (max > bmax)
    {
        cout << i << "max : " << max << endl;
        toStringInt(listeCodeBGC[i]);
        bmax = max;
    };

    // Update global minimum balance if this code is more balanced (shorter max run)
    if (max < bmin)
    {
        cout << i << "min : " << max << endl;
        toStringInt(listeCodeBGC[i]);
        bmin = max;
    }
}

/**
 * BGC with Population-based Pruning.
 * This version tracks how many elements remain in each Hamming layer
 * to ensure no layer is orphaned during the search.
 */
void BGC(int x)
{
    // Retourne tous les code de Beckett Gray sur n bits, x l'entier du code de Beckett Gray courant
    nbAppel++;
    if (s >= n2)
    {
        if (s > smax)
            smax++;

        // Si on a généré le code de Beckett Gray en entier
        listeCodeBGC.push_back(bgc);
        // On stocke dans listeCodeBGC et dans un fichier texte
        PE.stop();
        std::ofstream tmp("time.txt", std::ios::app);
        tmp << "\n"
            << PE.seconds();
        tmp.close();
        std::ofstream otp("codes.txt", std::ios::app);
        otp << "\n";
        for (int i : bgc)
        {
            otp << i;
            otp << " ";
        }
        otp.close();
        return;
    }
    if (s >= smax)
    {
        // On regarde la plus grande suite d'entiers de bgc pour l'instant générée
        if (s > smax)
            smax++;

        PE.stop();
        std::ofstream tmp("time.txt", std::ios::app);
        tmp << "\n"
            << smax << "\n"
            << PE.seconds();
        tmp.close();
        // On stocke dans un fichier texte
        std::ofstream otp("codes.txt", std::ios::app);
        otp << "\n"
            << smax << "\n";
        for (int i : bgc)
        {
            otp << i;
            otp << " ";
        }
        otp.close();
    }
    if (s > 1 && s - oldS[0] > nbMax[n])
        return;

    // --- NEW: LAYER ORPHAN PRUNING ---
    // If we moved to a higher weight but left unvisited nodes in the previous layer,
    // we can never go back to finish them. This branch is a dead end.
    if (V[x] > 1 and nbEltEns[V[x] - 1] == 0)
        return;

    int cmpt = 0; // Nombre de sommets isolés
    int xn = 0;   // Un sommet que l'on obtient en flippant un bit de x
    int bit = 0;  // Le bit qui change entre x et xn

    for (int i = 0; i < n; i++)
    {
        int fxi = flip[x][i];

        if (nombreVoisins[fxi] == 1 && avail[fxi])
        {
            // Si on a un voisins de x de degré 1
            cmpt++;
            if (cmpt > 1)
                return;

            xn = fxi;
            bit = n - i - 1;
        }
    }
    if (cmpt == 1 && (bi[xn][bit] == 1 || bit == old[0]))
    {
        // Si on a un sommet isolé voisin de x
        int oldS0 = oldS[0];

        if (bi[xn][bit] == 1)
        {
            // Si on passe d'un 0 à un 1, on ajoute bit dans l'odre d'entrée des 1 dans le code
            old.push_back(bit);
            oldS.push_back(s);
        }
        else
        {
            // Si on passe d'un 1 à un 0, on retire le premier élément de old, le deuxième bit de old devient le 1 rentré il y a le plus de temps
            old.erase(old.begin());
            oldS.erase(oldS.begin());
        }
        bgc.push_back(xn); // On ajoute xn au code de Beckett Gray
        for (int j = 0; j < n; j++)
            nombreVoisins[flip[xn][j]]--;

        avail[xn] = false;
        nbEltEns[V[xn]]--;

        s++; // On augmente la taille de bgc de 1

        BGC(xn); // On peut explorer les voisins de xn

        // On annule tous les changements fait précédemment pour pouvoir explorer d'autres solutions.
        nbEltEns[V[xn]]++;
        s--;

        for (int j = 0; j < n; j++)
            nombreVoisins[flip[xn][j]]++;

        if (bi[xn][bit] == 1)
        {
            old.pop_back();
            oldS.pop_back();
        }
        else
        {
            old.insert(old.begin(), bit);
            oldS.insert(oldS.begin(), oldS0);
        }
        avail[xn] = true;
        bgc.pop_back();
        return;
    }
    if (cmpt == 0)
    {
        // Si on a pas de sommet isolé
        for (int i = 0; i < n; i++)
        {
            // On parcours tous les voisins de x
            xn = flip[x][i];
            bit = n - i - 1;

            if (avail[xn] && (bi[xn][bit] == 1 || bit == old[0]))
            {
                // Le fonctionnement à identique à ci-dessus pour le cas du voisin isolé
                int oldS0 = (oldS.empty() ? 0 : oldS[0]);

                if (bi[xn][bit] == 1)
                {
                    old.push_back(bit);
                    oldS.push_back(s);
                }
                else
                {
                    old.erase(old.begin());
                    oldS.erase(oldS.begin());
                }
                bgc.push_back(xn);
                for (int j = 0; j < n; j++)
                    nombreVoisins[flip[xn][j]]--;
                s++;

                avail[xn] = false;
                nbEltEns[V[xn]]--;

                BGC(xn);

                nbEltEns[V[xn]]++;
                s--;

                for (int j = 0; j < n; j++)
                    nombreVoisins[flip[xn][j]]++;

                if (bi[xn][bit] == 1)
                {
                    old.pop_back();
                    oldS.pop_back();
                }
                else
                {
                    old.insert(old.begin(), bit);
                    oldS.insert(oldS.begin(), oldS0);
                }
                avail[xn] = true;
                bgc.pop_back();
            }
        }
        return;
    }
    return;
}

int main()
{
    // --- 1. SEARCH PARAMETERS ---
    nbMax = {0, 1, 2, 8, 16, 7, 10, 15, 20}; // Max consecutive 1s allowed
    n = 5;                                   // Dimension
    nbAppel = 0;                             // Reset call counter
    n2 = pow(2, n);                          // 2^n nodes
    smax = 0;                                // Best depth tracker

    // --- 2. GLOBAL STATE INITIALIZATION ---
    bgc = {0};         // Path starts at node 0
    listeCodeBGC = {}; // Container for found Gray codes

    // Setup visitation map (node 0 is already visited)
    avail.push_back(false);
    for (int i = 0; i < n2 - 1; i++)
        avail.push_back(true);

    // Initial degree for each node (n connections in an n-cube)
    for (int i = 0; i < n2; i++)
        nombreVoisins.push_back(n);

    // --- 3. LOOKUP TABLES & HAMMING WEIGHTS ---
    bi = listeBinaire(n); // Bit-vector lookup
    flip = listeFlip(n);  // Adjacency lookup
    s = 1;                // Current path length

    // Initialize the population count per Hamming layer
    for (int i = 0; i < n + 1; i++)
        nbEltEns.push_back(0);

    // Calculate Hamming weights and fill nbEltEns accordingly
    V = listeDes1(n);

    // --- 4. BECKETT FIFO QUEUE SETUP ---
    old = {};  // Active bit indices
    oldS = {}; // Timestamps of activation
    bmax = 0;
    bmin = 100000;

    // --- 5. PERFORMANCE MEASUREMENT ---
    PE.start(); // Start CPU cycle & time counters
    BGC(0);     // Launch recursive search
    PE.stop();  // Stop timers

    // --- 6. OUTPUTS ---
    std::cout << " nbr cycles : " << PE.cycles() << std::endl;
    std::cout << " nbr secondes : " << PE.seconds() << std::endl;
    std::cout << " nbr millisecondes : " << PE.milliseconds() << std::endl;
    cout << "Total codes: " << listeCodeBGC.size() << endl;
    cout << "Nombres d'appels: " << nbAppel << endl;

    PE.clear(); // Reset Perf monitor
    return 0;
}