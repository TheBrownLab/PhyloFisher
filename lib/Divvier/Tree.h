
#ifndef __TREE_HEADER
#define __TREE_HEADER

#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <regex>
#include <algorithm>
#include <iomanip>
#include <assert.h>
#include <list>
#include <cfloat>
#include "Random.h"
#include "Sequence.h"	// Some basic functions for file handling

using namespace::std;

// Some stuff originally in tools.h
string EatWhiteSpace(string line);
#define my_min(a,b) ((a)<(b)?(a):(b))
#define my_max(a,b) ((a)>(b)?(a):(b))
#define GET_MEM(var,type,size) var = new type[size]; assert(var != NULL);
#define DEL_MEM(var) if(var != NULL) { delete [] var; } var = NULL;
#define IFOR(iter, list) for(iter = list.begin(); iter!= list.end(); iter++)
#define FOR(i,n) for(i=0; i<n; i++)                             // Forwards for-loop
ostream& Error(string str = " ", ostream  &os = cout);
string int_to_string(int num);
string double_to_string(double num);
bool FlipBool(bool V);
bool FlipBin(int i);

template <class TOutVec> ostream& operator<<(ostream & os, vector <TOutVec> Vec) { int i; FOR(i,(int)Vec.size()) { os << Vec[i] << " "; } return os; }
template <class TIsIn> bool IsIn(TIsIn Val, vector <TIsIn> Vec) {
        if(Vec.empty()) { return false; }
        if(find(Vec.begin(),Vec.end(),Val) == Vec.end()) { return false; }
        return true;
}
template <class TSumType> TSumType Sum(vector <TSumType> *Vec) {
        TSumType Sum = 0;
        int i;
        FOR(i,(int)Vec->size()) { Sum += Vec->at(i); }
        return Sum;
}
template <class TComp> bool Compare(vector <TComp> *List1, vector <TComp> *List2) {
        int i;
        bool same = true;
        if(List1->size() != List2->size()) { return false; }
        FOR(i,(int)List1->size()) { if(fabs((double) (List1->at(i) - List2->at(i))) > 1.0E-6) { same = false; break; } }
        return same;
}





#define TREE_START 1			// Number that input trees first species start with
#define SET_BRANCH 0.1			// Default branch length

enum ENodeType	{ branch, leaf, root};

struct SSplit{ int BrLabel; vector <int> Left, Right; bool rootLeft = false; bool rootRight = false; };

class CNode		{
public:
    // Member variables
    ////////////////////////
	int m_iNoLinks;
	int m_iInternalNodeNum;
    vector <int> m_viLink;
    vector <int> m_viBranch;
    ENodeType m_NodeType;
    // Member functions
    //////////////////////
	// Note parent node ALWAYS should be specified in the last branchx and linkx where applicable
	// General constructor
	CNode(int NoLinks = -1,int *LinkList = NULL);
    // Constructor for internal nodes
    CNode(int linka, int linkb, int linkc, int brancha, int branchb, int branchc, int IntVal = -1);
    // Constructor for external nodes
    CNode(int linka, int brancha, int IntVal = -1);
	// Other constructor
	CNode(vector <int> Links, vector <int> Branches, int IntVal);
	// Copy constructor
	CNode(const CNode &Node);
    // Destructor function
    virtual ~CNode();
	// Member functions specific for bifurcating trees where linkc is defined as a parent
	void SetNode(int linka, int linkb, int linkc, int brancha, int branchb, int branchc, int IntVal);
	void SetNode(int la, int lb, int ba, int bb, int IntVal);
    void SetNode(int linka, int brancha, int IntVal);
	// General functions for assigning links (No parent assumed)
	void SetNode(int NoLinks, int *LinkList);
	// Cleaning operation
	void CleanNode();
	// Operator functions
    CNode &operator=(const CNode &);
};

ostream& operator<<(ostream& os, const CNode &Node);

class CTree
{
public:
    // Member functions
    /////////////////////////////////

    // Constructor functions
    CTree(string TREE, vector <string> Names, bool AllowFail = false, bool AllowSubTree = false);			// Basic constructor
    CTree();
	void CreateTree(string TREE,vector <string> Names,bool CheckVar = true, bool AllowFail = false,bool AllowSubTree = false);// Underlying construction function
	// Copy Constructor
	CTree(const CTree &Tree);
    // Destructor function
    virtual ~CTree();
	// Memory functions
	void CleanTree();
    void GetMemory(int NoSeq);

	// Access functions
	////////////////////////////////////////
	// Main tree
    inline bool IsRooted()  { return m_bRooted; }		// Returns whether the tree is rooted
    inline int Root()		{ return m_iRootNode; }			// Returns the root node index
	inline int NoNode()		{ return m_iNoNode; }		// Number of nodes in tree
	inline int NoBra()		{ return m_iNoBra; }		// Returns number of branches in tree
	inline int NoOptBra()	{ return m_iNoOptBra; }		// Returns the number of optimised branches in tree
	inline int NoSeq()		{ return m_iNoSeq; }		// Returns number of sequences in tree
	inline int StartCalc()	{ return m_iStartCalc; }	// Returns where the calculation starts
	bool OldTree()			{ return m_bOldTree; }		// Returns whether the tree is old or not
	bool SetOldTree(bool V)	{ m_bOldTree = V; return V; }	// Sets the m_bOldTree value
	int BestStartCalc();								// Sets the m_iStartCalc to the best value
	int SetStartCalc(int S) { m_iStartCalc=S;return S;}	// Forces m_iStartCalc to a value -- power user function because can cause problems
	bool FastCalcOK()		{ return m_bFastCalcOK; }	// Returns whether the fast calc is okay or needs resetting
	void SetFastCalcOK(bool V)	{ m_bFastCalcOK = V; }	// Sets the m_bFastCalcOK flag
	int SetOutBra(int V)	{ m_iOutBra = V; return V; }// Sets m_iOutBra
	void OutBra()			{ m_iOutBra = 1; }			// Sets m_iOutBra to output branches
	void NoOutBra()			{ m_iOutBra = 0; }			// Sets m_iOutBra to not output branches
	void OutBraNum()		{ m_iOutBra = 2; }			// Sets m_iOutBra to output branch numbers
	void OutName()			{ m_bOutName = true; }		// Set for outputting names
	void NoOutName()		{ m_bOutName = false; }		// Set for not outputting names
	void OutLabel()			{ m_bOutLabel = true; }		// Set for outputting labels on tree
	void NoOutLabel()		{ m_bOutLabel = false; }	// Set for not outputting labels on tree
	void CreateBranchLabels(bool Force = true);			// Create labels on branches that correspond with their numbers
	bool Valid() { return m_bValid; }					// Returns whether a valid constructor has been run
	void ForceReady() { m_bReady = true; }				// Makes the tree ready.
	vector <string> Names() { return m_vsName; }		// Returns the vector of names
	void SetNames(vector <string > NewNames, bool Overwrite = false);	// Set the names in the tree (e.g. for output)
	bool BranchLabels() { if(m_viBranchLabels.empty()) { return false; } return true; }
	vector <int> Labels() { return m_viBranchLabels; }
	int NoLabels() { return m_iNoBraLabels; }
	// Branch
	vector <double> Branches();					// Returns a vector of the branch lengths
	double B(int Branch);						// Returns the value of branch length
	double *OptimiserB(int Branch);				// Returns the pointer to the value to be used in optimisation
	double TreeLength();						// Returns the tree length;
	double SetB(int Branch, double Value, bool Update = false, bool Rescale = false);		// Set branch length Branch to Value
	double MulB(int Branch, double Value, bool Update = false, bool Rescale = false);		// Multiply branch length Branch by Value
	double AddB(int Branch, double Value, bool Update = false, bool Rescale = false);		// Add Value to branch length Branch
	double QuadB(int Branch);																// Get the average value of Branch and the (upto) 4 surrounding branches
	inline int BraLink(int B, int L) { return m_vBraLinks[B][L]; }						// Get the L link for branch B
	inline void ReplaceBraLink(int B,int L,int Val) { m_vBraLinks[B][L]=Val;};			// Set the L link for branch B to Val; returns Val;

	// Node
	inline int NoLinks(int N)	{ return (int) m_Node[N]->m_viLink.size(); }	// Returns the # links in a node
	int NodeLink(int N,int L)	{ return m_Node[N]->m_viLink[L]; }		// Returns link L of node N
	bool IsNodeLink(int N, int Val) { return IsIn(Val,m_Node[N]->m_viLink); }	// Returns whether a Val is in Node N m_viLink
	void ReplaceNodeLink(int N, vector<int>L);							// Replace m_viLink of node N with L;
	void ReplaceNodeLinkElement(int N, int Element, int Val);			// Replace a single element of m_viLink in Node N with Val
	int NodeBra(int N, int B)	{ return m_Node[N]->m_viBranch[B]; }	// Returns Branchlink B of Node N
	void ReplaceNodeBra(int N, vector <int>B);							// Replace m_viBranch of Node N with B
	void ReplaceNodeBraElement(int N, int Element, int Val);			// Replace a single element of m_viBranch in Node N with Val
	ENodeType NodeType(int N)	{ return m_Node[N]->m_NodeType; }		// Returns the NodeType of Node N
	bool NodeNull(int N) { if(m_Node[N] == NULL) { return true; } return false; }	// Returns whether a node is NULL or not
	int NoLeafLink(int N);												// Returns the number of leaf links for Node N
	void AssignNodeType(int N, ENodeType Type);							// Set the NodeType

	// Output functions
	friend ostream& operator<<(ostream& os, CTree &Tree);		// Standard output routine
	bool OutDetail(ostream &os = cout, bool ForceExit = false);	// Detailed output routine (nodes, branches, and tree)

	// Tree modification routines
	void MidpointRoot();											// Does midpoint rooting
	void AddRoot(int NodeLeft, double BrLeft, int NodeRight);		// Adds a node between NodeLeft and NodeRight
	void Unroot();											// If rooted, this function permanently unroots it
	void OrderNode(int NodeNum = -1,bool DoBraToo = true);	// Orders tree nodes to ascending numerical value (-1 does all)
	int CutBranch(int Branch,int Link);					// Cuts a branch: returns the branch on priority subtree where they used to be attached -- Subtrees are maintained; Link specifies which one is given priority
	int RemoveBraNode(int Node);						// Removes a node from a tree: returns the branch where link used to be attached -- NB: One Link must be set to -1
	int RemoveLeafNode(int RemNode);
	int AddSeq(int SeqNum, int Branch, double BranchProp =0.5);	// Adds a sequence to the tree
							// SeqNum = sequence added; Branch = branch added to; BranchProp = relative position in branch (from low to high node nums)
	int GetStart(bool replace = true);			// Find branch in tree from which to recurse
	double GetTreeLength(bool first = true, int NTo = -1, int NFr = -1);	// Get remaining stuff in a tree
	void BuildOriSubTree(CTree *T, vector <bool> NodesBool);										// Returns subtree from an array of bools describing which nodes it covers
	void BuildOriSubTree(CTree *T, vector <int> LeafMap, vector <int> NCover, vector <int> NFrom);	// Returns subtree from LeafMap and NodesCovered
	void BuildOriSubTree(CTree *T, vector <int> NodesCovered);										// Returns subtree from the nodes it covers
	void BuildBraLinks(bool Verify = true);					// Builds the branch links from nodes

	// Function to create a consistent output for comparison
	vector <int> ConstOut();
	int NodeDist(int Node1, int Node2, int NodeFrom = -1);		// Count number of branches separating nodes i and j;
	int BranchDist(int Br1, int Br2, bool AllowZero = false);	// Count the number of branches separating branches Br1 and Br2; !AllowZero means that very short branches will not be considered real branches
	// Functions for comparing trees using Robinson-Foulds distance
	int GetRFDist(CTree &Tree);			// Standard comparison between two trees of same number of taxa
	bool IsCompatible(CTree &SubTree);	// Compares tree <SubTree> with #Seq <= *this->#Seq to check whether they are compatible (High level function accepting tree objects)

	// Functions for adding sequences to a subtree tree based on an existing full tree (Greedy algorithm for maximising tree length)


	// Functions for testing properties of trees
	bool IsNode(int Node);  // Used for assert statements
	bool IsBra(int Branch); // Used for assert statements
	bool IsCutTree();		// Whether tree has had sequences removed...
	bool GoodBra(int Branch);	// Is an active branch in the tree
	bool GoodNode(int Node);	// Is an active node in the tree

	// Tree split-based functions
	vector <SSplit> BuildSplits();					// Build the set of splits associated with a tree. Current implementation always forces the rebuild
	SSplit GetSplit(int Bra, bool forceRebuild = false);		// Return the split set for branch Bra
	int BranchSets(int BranchNum, vector <int> &Left, vector <int> &Right);	// Find the sets of leaf sequences that the branch splits
	int FindBra(int Node_i, int Node_j);	// Find branch linking Nodes i and j, returns -1 if none
		// Returns the value of total number in the Left set (i.e. Left = m_ariBraLinks[0] )
	void GetBraSets(int NTo, int NFr, vector <int> &List, bool First = true);
	void OutSplits(ostream &os = cout);
	// Functions to get pairwise distances from a tree
	vector <double> GetTreePW();
	vector <double> GetAllTreePW();
	void PWDistSub(int NodeTo, int NodeFrom,vector <double> *d,bool DoInternalNodes = false);
	vector <double> GetPartialTreeDist(vector <int> LeafMap, vector <vector <int> > NBelow);
	vector <double> GetSubTreePW(vector <int> LeafMap, vector <vector <int> > NBelow, vector <double> BaseDist);
	vector <int> GetBranchPath(int x, int y, vector <int> current = vector<int>(), bool First = true);			// The the ordered list of branches from node y to x
	vector <int> GetNodePath(int x, int y, vector <int> current = vector<int>(), bool First = true);			// The the ordered list of branches from node y to x
	// Function to get all of the nodes of a certain depth from a particular node; bool GetLess == true, will also return leaf nodes when they fall within this range
	vector <int> GetNodesOfDepth(int InitNode, int NodeDepth, bool GetLess, vector <int> *NodesFrom = NULL, vector <int> *NodeCov = NULL, vector <double> *ExtBra = NULL ,int NodeFr = -1, bool First = true);
	// Centering point functions used for snap
	vector <int> BranchCP(int CP, int depth, vector <int> *NodeFr, vector <int> *NodesCovered = NULL, vector <double> *ExtBra = NULL);
	vector <int> NodeCP(int Node, int depth, vector <int> *NodeFr, vector <int> *NodesCovered = NULL, vector <double> *ExtBra = NULL);
	void ReplaceTreeCP(CTree *NT,vector <int> LeafMap,vector <int> NCover,bool VerifyBranchLinks = true);		// Overwrites the current tree (from CPs) and puts NewTree in its place
    // Operator= function
    void operator=(const CTree &);

	vector <vector <int> > GetKnotClusters(vector <bool> IntNodeChanges, int ChangeRad = 2);		// Function that takes an array[NoNodes()] of which branches change (or zero lengths, or whatever)
					// and returns a set of clusters of overlapping changes of SNAP radius = ChangeRad
private:
    // Private member variable
    ///////////////////////
	int m_iRootNode;			// The root node (-1 if not rooted)
	bool m_bRooted;				// Whether the tree is rooted
    int m_iNoNode;				// The # of nodes in tree
    int m_iNoBra;				// The # of branches in tree
	int m_iNoOptBra;			// The # of branches optimised in the current tree
    int m_iNoSeq;				// The # of sequences in tree
	int m_iStartCalc;			// Leaf node from which to start calculations
	vector <double> m_vdBra;	// The Branch length parameters
	vector <CNode *> m_Node;	// The nodes in the tree
//    CNode ** m_Node;			// The nodes in the tree
    bool m_bReady;				// Whether tree is ready
    int m_iOutBra;				// Out branch [0=no branches,1=branches,2=branch numbers]
	int m_bOutName;				// Out Names
	bool m_bOutLabel;			// Whether to output tree labels or not
	vector <string> m_vsName;	// The names in the sequence
	vector <vector <int> > m_vBraLinks; // the nodes linked to each branch [branch_number][links]
	bool m_bOldTree;			// Whether the tree has already been through a round of SNAP
	bool m_bFastCalcOK;			// Whether the tree is okay for FastCalc computations
	bool m_bValid;				// Flag to identify whether a valid constructor has been run
	vector <int> m_viBranchLabels;	// Labels for different branch types (e.g. parameter per branch models)
	int m_iNoBraLabels;				// Number of unique branch labels
	vector <SSplit> m_vSplits;	// Vector containing the splits on a tree

	// Private member functions
    ////////////////////////////////////////////////////////////

	// Find closest: gets closest nodes in a string
    void find_closest(string *tree, int *c1, int *c2, int *p,int n_p, double *bra,int *label, int *IntVal);
	void GetBraDetails(string s, int *node, double *bra, int *label);
	ostream& OutNode(int FromNodeNum, int ToNode, ostream &os);
	ostream& OutBranch(int ToNode, int FromNode, ostream &os);
	double DoBranch(string *tree,int *pos, int *IntVal);
	void BuildBranches(int NT = -1, int NF = -1);	// Builds the branchs of a tree from nodes recursively
	bool ValidateTree(bool AllowExit = true);

	// Some Split functions
	SSplit CalculateSplit(int Bra);						// Return the split set for branch Bra

	////////////////////////////////////////////
	// Rearrangement functions
	void BuildBraSA(CTree *T, int Seq, vector <CTree *> *TList, int OriBr, int SPRID, int MinDist, int MaxDist);
	void Branch_SA(CTree *T, int Seq, vector <CTree*> *TList, int First, int NTo, int NFr, int OriBr, int SPRID, int MinDist, int MaxDist);
	void DoSA(CTree *T, int Seq, vector <CTree*> *TList, int First, int NTo, int NFr, int Br, bool IsExtBra, int OriBr, int SPRID, int MinDist,int MaxDist);
};

bool IsSameTree(CTree *T1, CTree *T2);
void ExpandNode(int Node, string *String, int Stringpos, CTree *TREE);
int IsTreeGap(char Check);
// Functions for finding the greedy subtree
CTree FindGreedySubTree(CTree *FullTree, int NoSeq); 		// Driver function: uses the greedy algorithm to identify the optimal subtree
double TravAddGreedy(CTree *CurT, int To, int Fr, int Seq2Add, vector <double> *PWDist, int *BestBra); // Function to traverse a tree and test which position maximises distance
double GetDist(int Add, int a, int b,vector <double> *PWDist);	// Check how much improvement a given sequence added to the tree would provide
void GreedySeq2Tree(int Bra,int Seq2Add, CTree *CurTree, vector <double> *PWDists);  // Adds the sequence to the tree

bool SplitsCompatible(vector <SSplit> Splits1, int S1_seq, vector <SSplit> Splits2, int S2_seq);	// Low level function just comparing a set of splits
bool CompareSplit(SSplit S1, SSplit S2);															// Simple function for comparing splits

vector <string> ReadTreeNames(string Tree);

#endif
