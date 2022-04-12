#include <cmath>

#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"

#include "family.h"
#include "gene_family.h"
#include "WGDgc.h"

extern "C" {
#include "cafe.h"
}

using namespace std;

TEST_GROUP(WGDgc)
{
};

static pCafeTree create_tree(string tree)
{
    family_size_range range;
    range.min = range.root_min = 0;
    range.max = range.root_max = 15;

    char ctree[1000];
    strcpy(ctree, tree.c_str());
    return cafe_tree_new(ctree, &range, 0.01, 0);
}

TEST(WGDgc, binomialCoeff)
{
    LONGS_EQUAL(10, wgd::binomialCoeff(5, 2));
    LONGS_EQUAL(70, wgd::binomialCoeff(8, 4));
}

TEST(WGDgc, calcProbSamePara)
{
    DOUBLES_EQUAL(-0.224551, wgd::calcProbSamePara(log(0.05), 100, 8, 9), 0.0001);
}

TEST(WGDgc, calcProb)
{
    vector<double> para{ log(0.02), log(0.03) };
    DOUBLES_EQUAL(0.0160687, wgd::calcProb(para, 100, 8, 9), 0.0001);
}

TEST(WGDgc, calcProbOneBrach)
{
    vector<double> para{ log(0.02), log(0.03) };
    std::vector<pCafeNode> branchNode;
    std::vector<double> logLamlogMu;
    int nPos;
    int nFamily;
    int nLeaf;
    std::map<std::string, std::vector<int>> geneCountData;
    std::map<pCafeNode, wgd::Rate> wgdTab;
    std::map<pCafeNode, wgd::phylo_data> phyloMat;
    std::map<pCafeNode, wgd::Doom> MatDoomed;
    bool doomLeft;
    std::vector<std::vector<std::vector<double>>> Mat;

    wgd::calcProbOneBranch(branchNode, logLamlogMu, nPos, nFamily, nLeaf, geneCountData, wgdTab, phyloMat, MatDoomed, doomLeft, Mat);

    // DOUBLES_EQUAL(0.0160687, , 0.0001);
}

TEST(WGDgc, getPt_with_wgd)
{
    vector<int> leafcounts;
    vector<double> para{ log(0.02), log(0.03) };
    auto m = wgd::getPt(para, 1000, 5, -1, false, leafcounts, true, 0.01);
    LONGS_EQUAL(5, m.size());
    LONGS_EQUAL(5, m[1].size());
    DOUBLES_EQUAL(0.0001, m[3][3], 0.000001);
}

TEST(WGDgc, getPt_with_child)
{
    vector<int> leafcounts{3, 6, 4, 5, 5};
    vector<double> para{ log(0.02), log(0.03) };
    auto m = wgd::getPt(para, 1000, 7, 5, true, leafcounts, false, 0.01);
    LONGS_EQUAL(7, m.size());
    LONGS_EQUAL(5, m[1].size());
    DOUBLES_EQUAL(-0.526683, m[3][3], 0.000001);
}

TEST(WGDgc, getPt_with_internal)
{
    vector<int> leafcounts{ 3, 6, 4, 5, 5 };
    vector<double> para{ log(0.02), log(0.03) };
    auto m = wgd::getPt(para, 1000, 7, 5, false, leafcounts, false, 0.01);
    LONGS_EQUAL(7, m.size());
    LONGS_EQUAL(7, m[1].size());
    DOUBLES_EQUAL(-0.444404, m[3][3], 0.000001);
}

TEST(WGDgc, a)
{
    vector<string> species = { "A", "B", "C", "D" };
    string treestring = "(D:18.03,(C:12.06,(B:7.06,A:7.06):2.49):5.97);";
    //tre.string = "(D:{0,18.03},(C:{0,12.06},(B:{0,7.06},
    //    A:{0, 7.06}):{0, 2.49:wgd, 0 : 0, 2.50}):{0, 5.97});"

    auto phylo4d = create_tree(treestring);

    auto pfamily = cafe_family_init(species);
    gene_family gf[4];
    gf[0].id = "A";    gf[0].desc = "description";  gf[0].values = { 2, 2, 3, 1 };
    gf[1].id = "B";    gf[1].desc = "description";  gf[1].values = { 3, 0, 2, 1 };
    gf[2].id = "C";    gf[2].desc = "description";  gf[2].values = { 1, 0, 2, 2 };
    gf[3].id = "D";    gf[3].desc = "description";  gf[3].values = { 2, 1, 1, 1 };
    for (int i = 0; i<4; ++i)
        cafe_family_add_item(pfamily, gf[i]);

    auto a = wgd::processInput(phylo4d, 0.9);
    auto para = vector<double>{ log(0.01), log(0.02) };

    getLikGeneCount(para, a, pfamily, 8, 1 / 1.5, wgd::oneOrMore);

    cafe_tree_free(phylo4d);
    cafe_family_free(pfamily);
}