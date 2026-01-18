/*TER : Codes de Gray*/
/*Matté et Lucas*/

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <bitset>
#include "EvalPerf.hpp"

using namespace std;

int n;
int s;
int n_squared;
int smax;

vector<bool> avail;
vector<vector<int>> listeCodeBGC;
vector<int> p2;
vector<vector<int>> bi;
vector<vector<int>> flip;
vector<bool> l11;
vector<int> old;
vector<int> bgc;
vector<int> listev;

EvalPerf PE;
std::ofstream otp("codes.txt");

vector<int> bin_rep(int x, int n) {
    string b = bitset<32>(static_cast<unsigned long long>(x)).to_string(); 
    vector<int> res;
    
    for (int i = b.size() - n; i < b.size(); ++i)
        res.push_back(b[i] - '0');
    
    return res;
}

bool is_pw2(int n) {
    return n != 0 && (n & (n - 1)) == 0;
}

vector<int> pw2_list(int n) {
    vector<int> res;
    
    for (int i = 0; i < pow(2, n); ++i)
        if (is_pw2(i))
            res.push_back(i);
            
    return res;
}

vector<vector<int>> bin_list(int n) {
    vector<vector<int>> res;
    
    for (int i = 0; i < pow(2, n); ++i)
        res.push_back(bin_rep(i, n));
        
    return res;
}

vector<vector<int>> flip_list(int n) {
    vector<vector<int>> res;
    
    for (int i = 0; i < pow(2, n); ++i) {
        res.push_back(vector<int>());
        
        for (int j = 0; j < n; ++j)
            res[i].push_back(i ^ (1 << j));
    }
    
    return res;
}

vector<bool> listeUn1(int n) {
    vector<bool> res;
    
    int count;
    int num;
    
    for (int i = 0; i < pow(2, n); ++i) {
        count = 0;
        num = i;
        
        while (num > 0) {
            count += num & 1;
            num >>= 1;
        }
        
        if (count==1)
            res.push_back(true);
        else
            res.push_back(false);
    }
    return res;
}

int testSym(vector<int>& bgc){
    for (int i=0; i<listeCodeBGC.size(); i++){
        bool flag=true;
        for (int j=0; j<n2-1; j++)
            if (bgc[j+1]!=listeCodeBGC[i][n2-j-1]){
                flag=false;
                break;
            }
        if (flag)
            return 0;
    }
    return 1;
}

int plusGrandeDiff(){
    int cmpt=0;
    int cmpt1=0;
    int cmpt2=0;
    
    for (int i=0; i<listeCodeBGC.size(); i++){
        cmpt1=0;
        for (int j=0; j<n2-1; j++){
            if (abs(listeCodeBGC[i][j]-listeCodeBGC[i][j+1])>cmpt){
                cmpt=abs(listeCodeBGC[i][j]-listeCodeBGC[i][j+1]);
                cmpt1=0;
            } 
            else
                if (abs(listeCodeBGC[i][j]-listeCodeBGC[i][j+1])==cmpt)
                    cmpt1++;
        }
        if (cmpt1>cmpt2) 
            cmpt2=cmpt1;
    }
    cout<<cmpt2<<endl;
    return cmpt;
}

void toString(vector<bool> bonjour){
    for (bool b: bonjour)
        cout<< b << " ";
        
    cout<<endl;
}

void BGCHeuristique1(int x) {
    //lnode[s]+=1;
    /*if (s>smax){
        smax++;
        cout<<smax<<endl;
        otp << "\n";
        for (int i : bgc) {
            otp << i;
            otp << " "; // write data to the file
        }
    }*/
    
    if (s >= n_squared) {
        listeCodeBGC.push_back(bgc);
        //cout << "Plus grande différence: " << plusGrandeDiff() << endl;
        //cout<<testSym(bgc);
        otp << "\n";
        for (int i : bgc) {
            otp << i;
            otp << " "; // write data to the file
        }
        return;
    }
    
    if (l11[x]){
        bool flag = false;
        for (int i : p2)
            if (avail[i]) {
                flag = true;
                break;
            }
        if (!flag)
            return;
    }
    
    int cmpt = 0;
    int xn = 0;
    int bit = 0;
    
    for (int i = 0; i < n; i++) {
        int fxi=flip[x][i];
        
        if (listev[fxi] == 1 && avail[fxi]) {
            cmpt++;
            if (cmpt > 1)
                return;
            xn = fxi;
            bit = n-i-1;
        }
    }
    
    if (cmpt == 1 && (bi[xn][bit] == 1 || bit == old[0])) {
        if (bi[xn][bit] == 1)
            old.push_back(bit);
        else{
            old.erase(old.begin());
        }
        avail[xn] = false;
        bgc.push_back(xn);
        for (int j = 0; j < n; j++) {
            listev[flip[xn][j]]--;
        }
        s++;
        BGCHeuristique1( xn);
        s--;
        for (int j = 0; j < n; j++) {
            listev[flip[xn][j]]++;
        }
        if (bi[xn][bit] == 1)
            old.pop_back();
        else{
            old.insert(old.begin(),bit);
        }
        avail[xn] = true;
        bgc.pop_back();
        return;
    }
    
    if (cmpt == 0) {
        for (int i = 0; i < n; i++) {
            xn = flip[x][i];
            bit = n-i-1;
            
            if (avail[xn] && (bi[xn][bit] == 1 || bit == old[0])) {
                if (bi[xn][bit] == 1)
                    old.push_back(bit);
                else
                    old.erase(old.begin());
                     
                avail[xn] = false;
                bgc.push_back(xn);
                
                for (int j = 0; j < n; j++)
                    listev[flip[xn][j]]--;
                
                s++;
                BGCHeuristique1(xn);
                s--;
                
                for (int j = 0; j < n; j++)
                    listev[flip[xn][j]]++;
                
                if (bi[xn][bit] == 1)
                    old.pop_back();
                else
                    old.insert(old.begin(),bit);
                    
                avail[xn] = true;
                bgc.pop_back();
            }
        }
        return;
    }
return;
}



int main() {
    n = 5;
    n_squared = pow(2, n);
    smax=0;
    bgc = {0};
    listeCodeBGC = {};
    
    avail.push_back(false);
    for (int i=0; i<n2-1; i++)
    	avail.push_back(true);
    	
    for (int i=0; i<n2; i++)
        listev.push_back(n);
        
    p_squared = pw2_list(n);
    bi = bin_list(n);
    flip = flip_list(n);
    l11= listeUn1(n);
    s=1;
    
    int lnode[n_squared];
    for (int i=0; i<n_squared; i++)
        lnode[i]=0;
    
    old={};
    
    PE.start();
    BGCHeuristique1(0);
    PE.stop();
    
    otp.close();
    cout << " nbr cycles : " << PE.cycles() << endl;
    cout << " nbr secondes : " << PE.seconds() << endl;
    cout << " nbr millisecondes : " << PE.milliseconds() << endl;
    cout << "Total codes: " << listeCodeBGC.size() << endl;
    ofstream outputFile("output.txt"); // create a new output file or overwrite an existing one
    if (outputFile.is_open()) { // check if the file was opened successfully
        for (int i=0; i<n2; i++){
            outputFile << lnode[i];
            outputFile << ","; // write data to the file
        }
        
        outputFile.close(); // close the file when done
        cout << "Data was written to output.txt\n";
    }
    else
        cerr << "Error opening file\n";
    PE.clear();
    cout << "Plus grande différence: " << plusGrandeDiff() << endl;
    return 0;
}
