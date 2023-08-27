#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "GLPK/include/glpk.h"

#include "LEDA/core/string.h"
#include "LEDA/core/tuple.h"

using namespace std;

struct STG
{
    int numInputs, numOutputs, numEdges, numStates;
    map<string, int> states;
    vector<vector<int>> edges;

    int & edge(string source, string target) { return edges[states[source]][states[target]]; }
    int & outEdge(string source) { return edges[states[source]][numStates + 1]; }
};

void readFile(const string inputFile, STG & stg)
{
    fstream input (inputFile, ios::in);
    for(int i = 0; i < 6; i++)
    {
        string tmp;
        input >> tmp;
        if(tmp == ".i")
            input >> stg.numInputs;
        else if(tmp == ".o")
            input >> stg.numOutputs;
        else if(tmp == ".p")
            input >> stg.numEdges;
        else if(tmp == ".s")
            input >> stg.numStates;
        else
            getline(input, tmp);
    }

    stg.edges.resize(stg.numStates + 1, vector<int> (stg.numStates + 2, 0));

    for(int i = 0; i < stg.numEdges; i++)
    {
        string in, source, target, out;
        input >> in >> source >> target >> out;
        if(stg.states.find(source) == stg.states.end())
            stg.states[source] = stg.states.size() + 1;
        if(stg.states.find(target) == stg.states.end())
            stg.states[target] = stg.states.size() + 1;
        int num = 1;
        for(const char & bit : in)
        {
            if(bit == '-')
                num *= 2;
        }
        stg.edge(source, target) += num;
        stg.outEdge(source) += num;
    }
    input.close();

    // cout << "Testing in readFile():\n"; cout.flush();
    // while(true)
    // {
    //     string source, target;
    //     cin >> source;
    //     if(source == "-1")
    //         break;
    //     cin >> target;
    //     cout << stg.edge(source, target);
    // }

    // for(const auto & source : stg.states)
    // {
    //     cout << source.first << " ";
    //     for(const auto & target : stg.states)
    //     {
    //         cout << stg.edge(source.first, target.first) << " ";
    //     }
    //     cout << stg.outEdge(source.first);
    //     cout << endl;
    // }
}

void probabilityCalculation(STG stg, vector<vector<double>> & weight)
{
    glp_prob *lp = glp_create_prob();
    glp_add_rows(lp, stg.numStates + 1);
    glp_add_cols(lp, stg.numStates);
    for(int i = 1; i <= stg.numStates; i++)
    {
        glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
        glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
    }
    glp_set_row_bnds(lp, stg.numStates + 1, GLP_FX, 1.0, 1.0);

    int numConstraints = 0;
    int ia[(stg.numStates + 1) * stg.numStates + 1], ja[(stg.numStates + 1) * stg.numStates + 1];
    double ar[(stg.numStates + 1) * stg.numStates + 1];
    for(const auto target : stg.states)
    {
        for(const auto source : stg.states)
        {
            if(target == source)
            {
                numConstraints++;
                ia[numConstraints] = target.second;
                ja[numConstraints] = source.second;
                ar[numConstraints] = -1.0 + double(stg.edge(source.first, target.first)) / stg.outEdge(source.first);
                // cout << "A\n";
                // cout << numConstraints << " | ia: " << target.first << " ja: " << source.first << " ar: " << -1.0 + double(stg.edge(source.first, target.first)) / stg.outEdge(source.first) << endl; cout.flush();
                // cout << numConstraints << " | ia: " << ia[numConstraints] << " ja: " << ja[numConstraints] << " ar: " << ar[numConstraints] << endl; cout.flush();
            }
            else if(stg.edge(source.first, target.first) > 0)
            {
                numConstraints++;
                ia[numConstraints] = target.second;
                ja[numConstraints] = source.second;
                ar[numConstraints] = double(stg.edge(source.first, target.first)) / stg.outEdge(source.first);
                // cout << "B\n";
                // cout << numConstraints << " | ia: " << target.first << " ja: " << source.first << " ar: " << double(stg.edge(source.first, target.first)) / stg.outEdge(source.first) << endl; cout.flush();
                // cout << numConstraints << " | ia: " << ia[numConstraints] << " ja: " << ja[numConstraints] << " ar: " << ar[numConstraints] << endl; cout.flush();
            }
        }
        numConstraints++;
        ia[numConstraints] = stg.numStates + 1;
        ja[numConstraints] = target.second;
        ar[numConstraints] = 1.0;
        // cout << "C\n";
        // cout << numConstraints << " | ia: " << stg.numStates + 1 << " ja: " << target.first << endl; cout.flush();
        // cout << numConstraints << " | ia: " << ia[numConstraints] << " ja: " << ja[numConstraints] << " ar: " << ar[numConstraints] << endl; cout.flush();
    }
    
    glp_load_matrix(lp, numConstraints, ia, ja, ar);
    glp_exact(lp, NULL);
    double solution[stg.numStates + 1];
    double min = 1;
    for(const auto & state : stg.states)
    {
        solution[state.second] = glp_get_col_prim(lp, state.second);
        min = (min < solution[state.second]) ? min : solution[state.second];
    }
    glp_delete_prob(lp);

    for(const auto & source : stg.states)
    {
        for(const auto & target : stg.states)
        {
            if(source == target)
                break;
            double sum = 0;
            if(stg.edge(source.first, target.first) > 0)
                sum += solution[source.second] * stg.edge(source.first, target.first) / stg.outEdge(source.first);
            if(stg.edge(target.first, source.first) > 0)
                sum += solution[target.second] * stg.edge(target.first, source.first) / stg.outEdge(target.first);
            if(sum > 0)
                weight[source.second][target.second] = weight[target.second][source.second] = sum * stg.numStates / sqrt(min);
        }
    }

    cout << "-----\n";
    for(const auto & state : stg.states)
        cout << state.first << ":\t" << solution[state.second] << endl;
    // cout << "Min: " << min << endl;
    cout << "-----\n";
    for(const auto & row : stg.states)
    {
        for(const auto & col : stg.states)
        {
            if(row == col)
                break;
            double prob = weight[row.second][col.second];
            if(prob > 0)
                cout << row.first << " <-> " << col.first << ":\t" << prob << endl;
        }
    }
}

int stateAssignment(STG stg, vector<vector<double>> & weight, vector<string> & stateCode)
{
    int length = 0;
    vector<unordered_set<int>> stateGroup;
    unordered_set<int> tmp;
    for(const auto & state : stg.states)
        tmp.emplace(state.second);
    stateGroup.emplace_back(tmp);
    tmp.clear();

    int round = 0;
    while(stateGroup.size() < stg.states.size())
    {
        cout << "   [   Round " << ++round << "   ]   \n";
        glp_prob *mlp = glp_create_prob();
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        glp_add_cols(mlp, stg.numStates);
        for(const auto & i : stg.states)
        {
            cout << "Set col: " << i.first << " -> " << i.second << endl;
            glp_set_col_kind(mlp, i.second, GLP_BV);
        }
        // for(int i = 1; i <= stg.numStates; i++)
        //     glp_set_col_kind(mlp, i, GLP_BV);
        
        int rowCount = 0, numConstraints = 0;
        int ia[stg.numStates * 2 + stg.numEdges * 12 + 2], ja[stg.numStates * 2 + stg.numEdges * 12 + 2];
        double ar[stg.numStates * 2 + stg.numEdges * 12 + 2];
        cout << "Size of Group: " << stateGroup.size() << endl;
        int numGroups = 0;
        for(const auto & group : stateGroup)
        {
            cout << "Group " << ++numGroups << " | Size " << group.size() << endl;
            if(group.size() > 1)
            {
                rowCount += 2;
                glp_add_rows(mlp, 2);
                int constA = min(int(ceil(log2(stg.numStates))) + 0, stg.numStates);
                if(constA == length + 1)
                {
                    cout << "Set row: " << rowCount - 1 << " lb: " << 1 << " ub: " << pow(2, constA - 1 - length) << endl;
                    glp_set_row_bnds(mlp, rowCount - 1, GLP_FX, 1.0, pow(2, constA - 1 - length));
                    cout << "Set row: " << rowCount << " lb: " << 1.0 - double(group.size()) << " ub: " << pow(2, constA - 1 - length) - group.size() << endl;
                    glp_set_row_bnds(mlp, rowCount, GLP_FX, 1.0 - double(group.size()), pow(2, constA - 1 - length) - group.size());
                }
                else
                {
                    cout << "Set row: " << rowCount - 1 << " lb: " << 1 << " ub: " << pow(2, constA - 1 - length) << endl;
                    glp_set_row_bnds(mlp, rowCount - 1, GLP_DB, 1.0, pow(2, constA - 1 - length));
                    cout << "Set row: " << rowCount << " lb: " << 1.0 - double(group.size()) << " ub: " << pow(2, constA - 1 - length) - group.size() << endl;
                    glp_set_row_bnds(mlp, rowCount, GLP_DB, 1.0 - double(group.size()), pow(2, constA - 1 - length) - group.size());
                }
                for(const auto & state : group)
                {
                    numConstraints++;
                    ia[numConstraints] = rowCount - 1;
                    ja[numConstraints] = state;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << endl;

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = state;
                    ar[numConstraints] = -1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << endl;
                }
            }
        }
        int colCount = stg.numStates;
        int numEdges = 0;
        for(const auto & source : stg.states)
        {
            for(const auto & target : stg.states)
            {
                if(source == target)
                    break;
                if(weight[source.second][target.second] > 0)
                {
                    numEdges++;
                    colCount++;
                    cout << "Edge found! " << source.first << ":" << source.second << " <-> " << target.first << ":" << target.second;
                    cout << "\tWeight: " << weight[source.second][target.second] << "\tSet col: " << colCount << endl;
                    glp_add_cols(mlp, 1);
                    glp_set_col_kind(mlp, colCount, GLP_BV);
                    // glp_set_col_bnds(mlp, colCount, GLP_DB, -1, 2);
                    glp_set_obj_coef(mlp, colCount, weight[source.second][target.second]);

                    rowCount++;
                    cout << "Set row: " << rowCount << " ub: " << 0 << endl;
                    glp_add_rows(mlp, 1);
                    glp_set_row_bnds(mlp, rowCount, GLP_UP, 0, 0.0);
                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = source.second;
                    ar[numConstraints] = -1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << " | ";

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = target.second;
                    ar[numConstraints] = -1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << " | ";

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = colCount;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << endl;

                    rowCount++;
                    cout << "Set row: " << rowCount << " ub: " << 2 << endl;
                    glp_add_rows(mlp, 1);
                    glp_set_row_bnds(mlp, rowCount, GLP_UP, 0, 2.0);
                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = source.second;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << " | ";

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = target.second;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << " | ";

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = colCount;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << endl;

                    rowCount++;
                    cout << "Set row: " << rowCount << " lb: " << 0 << endl;
                    glp_add_rows(mlp, 1);
                    glp_set_row_bnds(mlp, rowCount, GLP_LO, 0.0, 0);
                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = source.second;
                    ar[numConstraints] = -1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << " | ";

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = target.second;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << " | ";

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = colCount;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << endl;

                    rowCount++;
                    cout << "Set row: " << rowCount << " lb: " << 0 << endl;
                    glp_add_rows(mlp, 1);
                    glp_set_row_bnds(mlp, rowCount, GLP_LO, 0.0, 0);
                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = source.second;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << " | ";

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = target.second;
                    ar[numConstraints] = -1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << " | ";

                    numConstraints++;
                    ia[numConstraints] = rowCount;
                    ja[numConstraints] = colCount;
                    ar[numConstraints] = 1;
                    cout << numConstraints << ": " << ia[numConstraints] << " " << ja[numConstraints] << " " << ar[numConstraints] << endl;
                }
            }
        }
        rowCount++;
        glp_add_rows(mlp, 1);
        glp_set_row_bnds(mlp, rowCount, GLP_UP, 0, stg.numStates);
        for(const auto & state : stg.states)
        {
            numConstraints++;
            ia[numConstraints] = rowCount;
            ja[numConstraints] = state.second;
            ar[numConstraints] = 2;
        }
        glp_load_matrix(mlp, numConstraints, ia, ja, ar);
        glp_intopt(mlp, &parm);
        cout << "-----\n";
        cout << "Objective Function Value: " << glp_mip_obj_val(mlp) << endl;
        int num0 = 0, num1 = 0;
        for(const auto & state : stg.states)
        {
            cout << state.first << "\t" << state.second << "\t" << glp_mip_col_val(mlp, state.second) << endl;
            stateCode[state.second][length] = char(int(glp_mip_col_val(mlp, state.second)) + '0');
            if(int(glp_mip_col_val(mlp, state.second)) == 0)
                num0++;
            else
                num1++;
        }
        cout << "Number of 0 assigned: " << num0 << " | Number of 1 assigned: " << num1 << endl;
        for(int z = 1; z <= numEdges; z++)
        {
            cout << "XOR Constraint " << stg.numStates + z << ": " << glp_mip_col_val(mlp, stg.numStates + z) << endl;
        }

        glp_delete_prob(mlp);

        const int original = stateGroup.size();
        for(int i = 0; i < original; i++)
        {
            if(stateGroup[i].size() == 1)
                continue;
            unordered_set<int> & nowGroup = stateGroup[i];
            for(auto iter = nowGroup.begin(); iter != nowGroup.end(); )
            {
                if(stateCode[*iter][length] == '0')
                    iter++;
                else
                {
                    tmp.emplace(*iter);
                    iter = nowGroup.erase(iter);
                }
            }
            if(!tmp.empty())
            {
                stateGroup.emplace_back(tmp);
                tmp.clear();
            }
        }

        for(const auto & source : stg.states)
        {
            for(const auto & target : stg.states)
            {
                if(source == target)
                    break;
                if(stateCode[source.second][length] != stateCode[target.second][length])
                {
                    weight[source.second][target.second] *= 2;
                    weight[target.second][source.second] *= 2;
                }
            }
        }

        length++;
        cout << "-----\n";
        cout << "Size of State Group: " << stateGroup.size() << endl;
        for(const auto & group : stateGroup)
        {
            cout << " -";
            for(const auto & i : group)
                cout << " " << i;
            cout << endl;
        }
        cout << endl;
    }

    cout << "Final Results:\n";
    cout << "Used Length / Total Length:\t" << length << " / " << stateCode[0].size() << endl;
    for(const auto & state : stg.states)
        cout << state.first << "\t: " << state.second << "\t| " << stateCode[state.second] << endl;

    return length;
}

void writeFile(const string inputFile, const string outputFile, const STG stg, const vector<string> & stateCode, const int length)
{
    const string name = inputFile.substr(string("./benchmarks/").size(),
                                         inputFile.find(".kiss") - string("./benchmarks/").size());
    fstream input (inputFile, ios::in), output (outputFile, ios::out);
    output << ".model " << name << endl;
    string tmp;
    while(getline(input, tmp))
        output << tmp << endl;
    input.close();

    for(const auto & state : stg.states)
        output << ".code " << state.first << " " << stateCode[state.second].substr(0, length) << endl;
    output << ".end" << endl;
    output.close();
}

int main(int argc, char * argv[])
{
    string inputFile = argv[1];
    string outputFile = string("./results/") +
                        inputFile.substr(string("./benchmarks/").size(),
                                         inputFile.find(".kiss") - string("./benchmarks/").size()) +
                        string(".blif");
    
    cout << endl << string (75, '=') << endl;
    cout << "Reading File ...\n";
    STG stg;
    readFile(inputFile, stg);
    cout << endl << string (75, '=') << endl;

    cout << "Calculating State Probability ...\n";
    vector<vector<double>> weight (stg.numStates + 1, vector<double> (stg.numStates + 1, 0.0));
    probabilityCalculation(stg, weight);
    cout << endl << string (75, '=') << endl;

    cout << "State Assignment ...\n";
    vector<string> stateCode (stg.numStates + 1, string (stg.numStates, '-'));
    int length = stateAssignment(stg, weight, stateCode);
    cout << endl << string (75, '=') << endl;

    cout << "Writing File ...\n";
    writeFile(inputFile, outputFile, stg, stateCode, length);
    cout << endl << string (75, '=') << endl;
}