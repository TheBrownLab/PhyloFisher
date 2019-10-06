// Code for classes CTree and CNode
/////////////////////////////////////////////
// WARNING: This code has been heavily tampered with
// Some of the entry and exit condition of CTree function are *NOT* tested


// Include files
#include "Tree.h"
#define TREE_DEBUG 0


#if DO_MEMORY_CHECK
extern CMemChecker memory_check;
#endif

// Firstly describe CNode
/////////////////////////////////////////////////

// General constructor
////////////////////////////////////////////
CNode::CNode(int NoLinks,int *LinkList)	{
#if DO_MEMORY_CHECK
	memory_check.CountCNode++;
#endif
	SetNode(NoLinks,LinkList); }
CNode::CNode(int la, int lb, int lc, int ba, int bb, int bc,int IntVal) {
#if DO_MEMORY_CHECK
	memory_check.CountCNode++;
#endif
SetNode(la,lb,lc,ba,bb,bc,IntVal); }
CNode::CNode(int la, int lb, int IntVal) {
#if DO_MEMORY_CHECK
	memory_check.CountCNode++;
#endif
SetNode(la,lb,IntVal); }
CNode::CNode(vector <int> L, vector <int>B, int IntVal)	{
#if DO_MEMORY_CHECK
	memory_check.CountCNode++;
#endif
	assert(L.size() == B.size());
	switch(L.size())	{
	case 1: SetNode(L[0],B[0],IntVal); break;
	case 2: SetNode(L[0],L[1],B[0],B[1],IntVal); break;
	case 3: SetNode(L[0],L[1],L[2],B[0],B[1],B[2],IntVal); break;
	default: Error("Unknown type of node in constructor...");
	}
}
CNode::CNode(const CNode &Node)	{int i;
#if DO_MEMORY_CHECK
	memory_check.CountCNode++;
#endif
	CleanNode();
    m_NodeType = Node.m_NodeType;
	FOR(i,(int)Node.m_viBranch.size())	{ m_viBranch.push_back(Node.m_viBranch[i]); }
	FOR(i,(int)Node.m_viLink.size())		{ m_viLink.push_back(Node.m_viLink[i]); }
}
// Destructor function
///////////////////////////////////////////
CNode::~CNode() {
#if DO_MEMORY_CHECK
	memory_check.CountCNode--;
#endif
	CleanNode(); }

// SetNode functions for bifurcating trees
// These are appallingly verbose, but work...
///////////////////////////////////////////////
void CNode::SetNode(int la,int lb, int IntVal)	{	// Leaf node
	CleanNode();
	m_NodeType = leaf;
	m_viLink.push_back(la);
	m_viBranch.push_back(lb);
	m_iInternalNodeNum = IntVal;
}

// Internal node
void CNode::SetNode(int la, int lb, int lc, int ba, int bb, int bc, int IntVal)	{
	CleanNode();
    m_NodeType = branch;
	m_viLink.push_back(la); m_viLink.push_back(lb); m_viLink.push_back(lc);
	m_viBranch.push_back(ba); m_viBranch.push_back(bb); m_viBranch.push_back(bc);
	m_iInternalNodeNum = IntVal;
}

// Root node
void CNode::SetNode(int la, int lb, int ba, int bb, int IntVal) {
	CleanNode();
	m_NodeType = root;
	m_viLink.push_back(la); m_viLink.push_back(lb);
	m_viBranch.push_back(ba); m_viBranch.push_back(bb);
	m_iInternalNodeNum = IntVal;
}

/////////////// General links creation function ////////////////////
////////////////////////////////////////////////////////////////////

void CNode::SetNode(int NoLinks, int *LinkList)	{
	int i;
	CleanNode();
	// Specify leaf node
	if(NoLinks < 1) {
		m_NodeType = leaf;
	} else	if(NoLinks == 2) {
		FOR(i,NoLinks) { m_viLink.push_back(LinkList[i]); }
		m_NodeType = root;
	} else {	// Specify branch node
		FOR(i,NoLinks) { m_viLink.push_back(LinkList[i]); }
		m_NodeType = branch;
	}
}

// Clean Node
////////////////////////////////////////////////////////
void CNode::CleanNode()	{
	m_iInternalNodeNum = -1;
	m_viBranch.clear();
	m_viLink.clear();
}

// Stream operators for CNode
///////////////////////////////////////////
ostream& operator<<(ostream& os, const CNode &Node)	{
    int i;
	switch(Node.m_NodeType)	{
	case branch:	os << "\nInternal node "; break;
	case leaf:	os << "\nExternal node "; break;
	case root: os << "\nRoot node "; break;
	default:		os << "\n\nUnknown node type operator<<:(\n\n"; exit(-1);
	}
	if(Node.m_iInternalNodeNum >= 0) { os << "[ "<<Node.m_iInternalNodeNum<<" ] "; }
    os << " \tLinks:  "; FOR(i,(int)Node.m_viLink.size()) { os << Node.m_viLink[i] << "  "; }
    os << " \tBranches:  "; FOR(i,(int)Node.m_viBranch.size()) { os << Node.m_viBranch[i] << "  "; }
    return os;
}

// Node related operators
////////////////////////////////////////////////////////

CNode &CNode::operator=(const CNode & Node) {
    int i;
	CleanNode();
    m_NodeType = Node.m_NodeType;
	m_iInternalNodeNum = Node.m_iInternalNodeNum;
	FOR(i,(int)Node.m_viBranch.size())	{ m_viBranch.push_back(Node.m_viBranch[i]); }
	FOR(i,(int)Node.m_viLink.size())		{ m_viLink.push_back(Node.m_viLink[i]); }
	return *this;
}

// **********************************************************
// 	Description of CTree
// **********************************************************

// Blank constructor
///////////////////////////////////////////////////
CTree::CTree() {
#if DO_MEMORY_CHECK
        memory_check.CountCTree++;
#endif
    m_bReady = false ;m_iOutBra=0;
        m_iNoNode = m_iNoBra = m_iNoOptBra = m_iNoSeq = m_iStartCalc = 0;
        m_bOutName = false; m_bOldTree = false; m_bOutLabel = false; m_bRooted = false; m_iRootNode = -1;
}


// Top level constructors
CTree::CTree(string TREE, vector <string> Names, bool AllowFail, bool AllowSubTree)	{
#if DO_MEMORY_CHECK
	memory_check.CountCTree++;
#endif
	m_iRootNode = -1; m_bRooted = false; m_bValid = false; CreateTree(TREE,Names,true,AllowFail,AllowSubTree); }		// Basic constructor

int IsTreeGap(char Check)	{
	if(Check!='('&&Check!=')'&&Check!=';'&&Check!=','&&Check!=':')		{ return 1; }
	return 0;
}

// CTree copy constructor
/////////////////////////////////////////////////////
CTree::CTree(const CTree &Tree)	{
#if DO_MEMORY_CHECK
	memory_check.CountCTree++;
#endif
	m_bValid = false;
    m_bReady = false ;m_iOutBra=false; m_bOldTree = Tree.m_bOldTree;
	m_iNoNode = m_iNoBra = 0;
	*this = Tree;
}

// Underlying Constructor function
/////////////////////////////////////////////////////
// Takes two arguments: char tree[], and int NoSeq
// tree[] is the tree is bracket form eg. ((1,2)3,4)
// NOTE1:  The nnumber of brackets determines whether the tree is classified as rooted
//          or unrooted.  Brackets == NoSeq-1 => rooted, Brackets == NoSeq-2 => unrooted
// NOTE2:  The two species case is a special case.  When two sequences are observed the tree
//			is input in the form of a rooted tree.
// NOTE3:  !!!TODO:!!! This code expects reversible models to transform rooted tree to unrooted tree
// 			for efficient calculation
// NOTE4:  This routine is wasteful of space as it takes an extra node rooted trees regardless of their shape
// NoSeq is the number of leaf nodes in the tree
/////////////////////////////////////////////////////
void CTree::CreateTree(string Tree, vector <string> Names, bool CheckVar, bool AllowFail,bool AllowSubTree) {
    int i,j,pRight=0,pLeft=0,NextPointer=0,Parent, IntVal = -1, mem_seq, countBra, SubSeq;
    int NoSeq = Names.size();	// Get the number of sequences from the names provide
    double TempBranch[3];		// Stores branch lengths
	int TempLabel[3];			// Stores branch labels
    bool *CheckAlloc=NULL;		// Checks node allocation
	string TempTree;
	vector <string> Toks;
	list <int> temp_l;		// Space for reorganising branch labels
	vector <int> temp_l2;	// More space
	m_vsName = Names;
	// Do some checks on the string Tree
	// 1. Check ends in ';'
	if(Tree.find(";",(int)Tree.size()-1) == string::npos) { Tree += ";"; }
	// Count the '(' & ')' to calculate # of sequences and remove all spaces
	pLeft = pRight = j = 0;
	FOR(i,(int)Tree.size() - 1) {
		if(Tree[i] == '(') { pLeft++;} if(Tree[i] == ')') { pRight ++; }	// Count number of parentheses
		if((Tree[i] == '(' || Tree[i] == ',') && isalnum(Tree[i+1]) != 0) { j++; }	// Count hte number of species
	}

	////////////////////////////////////////////////////////////////
	// Do some error checking and check the number of species
	// 1. Left and Right parentheses
	if(pRight != pLeft) {
		cout << "\nUnequal number of pLeft(="<<pLeft<<") and pRight(="<<pRight<<") parentheses in tree:\n"<<Tree<<"\n\n";
		if(AllowFail) { return; } else { Error(); }
	}
	// 2. Number of species agree
	if(NoSeq > 0 && j != NoSeq && AllowSubTree == false) {
		if(CheckVar) { Error("\nCTree: Disagreement between number of species in file header ("+int_to_string(NoSeq)+") and number of species in file ("+int_to_string(j)+")"); }
		else {
			// Make sure difference between input tree and space reserved
			mem_seq = NoSeq; NoSeq = j;
		}
	} else {
		if(AllowSubTree) { assert(NoSeq > j); mem_seq = NoSeq; SubSeq = j; }
		else { mem_seq = NoSeq = SubSeq = j; }
	}
	// 3. Number of internal nodes is correct
	if(NoSeq == 2) { m_bRooted = false; }
	else if(pRight!=(NoSeq-1) && NoSeq>2) {		// For unrooted trees #internal nodes = #seq - 2
		m_bRooted = false;
	} else if(pRight!=(NoSeq-2) && NoSeq>2) {	// For rooted trees #internal nodes = #seq - 1
		m_bRooted = true;
	} else {
		cout << "\nIncorrect number of nodes for sequence number("<<NoSeq<<"):\n\n";
	    cout << pRight << " internal nodes, represented by '(', instead of ";
        cout << (NoSeq-2) << "\nTree:\n" << Tree;
		cout << "\n\n";
		if(AllowFail) { return; } Error();
	}
	// 4. Ends in ;
	if(Tree[(int) ( Tree.size() -1 )] != ';') { cout << "\nNo ; at end of tree:\n"<<Tree<<"\n\n"; if(AllowFail) { return; } Error();	 }
	///////////////////////////////////////////////////////////////////
	// Initialise
    m_bReady = false; m_iOutBra=m_bOutName=false; m_bOldTree = false; m_bOutLabel = false;
	GetMemory(mem_seq);
	m_viBranchLabels.assign(m_iNoBra,0);
	m_iStartCalc = 0;
	m_bFastCalcOK = false;
	CheckAlloc = new bool[m_iNoNode+1]; assert(CheckAlloc != NULL);
    for(i=0;i<m_iNoNode+1;i++) { CheckAlloc[i] = false; }
    // Get new list of names and their index
    vector<tuple<string ,int >> newNames;
    i = 0;
    for(auto name : Names) {
    	newNames.push_back(tuple<string,int>(name,i++));
    }
    // Sort names by length so never a problem with one being a substring of the other
    sort(newNames.begin(), newNames.end(),[](const auto a, const auto b) {
    	return get<0>(a).size() > get<0>(b).size();
    });
	// Put the numbers into the tree
    for(auto name : newNames) {
    	int pos = Tree.find(get<0>(name));
    	if(pos == string::npos) { cout << "\nERROR: Sequence name not found in tree: " << get<0>(name) << "\n"; exit(-1); }
    	Tree.replace(pos, get<0>(name).size(), int_to_string(get<1>(name) + TREE_START));
    }
	FOR(i,(int)Tree.size()) {
		// Check it's a proper branch length
		if(Tree[i] == ':') {
			i++;
			for(j = i;j<(int)Tree.size();j++) {
				if(Tree[j] == ',' || Tree[j] == ')') { i = j; break; }
			}
		}
		if(isalpha(Tree[i])) { Error("\nError in CTree::CreateTree(): tree has names not passed to the function...\nTree: " + Tree + "\n\n"); }
	}
	TempTree = Tree;

    if(SubSeq>2)	{		// Prepare first bit of trNextPointeree if more than two species
        // Loop to allocate nodes pLeft, pRight, Next are the old child0 child1 and parent nodes
		// these are allocated to m_Node->link[0], [1] and [2] respectively  This is used to some effect
		// when allocated branches and so on be careful not to screw it up
    	if(m_bRooted) { countBra = SubSeq-1; } else { countBra = SubSeq-3; }	// Note SubSeq used to allow subtrees
    	// Note the last node is done in this routine for rooted trees and finalised in the next section
        for(NextPointer = m_iNoSeq, i=0;i<countBra;i++)	{	// Get the non-trifurcating bits of tree
//        	cout << "\nProcessing: " << TempTree << flush;
			for(j=0;j<3;j++) { TempBranch[j] = SET_BRANCH; }
            // Find open bracket details and process
            find_closest(&TempTree,&pLeft,&pRight,&Parent,NextPointer, TempBranch, TempLabel,&IntVal);
			// Initialise the required nodes
			if(pLeft < m_iNoSeq) 	{	// If pLeft is leaf then initialise
				if(CheckAlloc[pLeft] != false)	{ cout << "\nNode " << pLeft << " multiply allocated.\nUsually caused by same sequence occuring multiple times in tree\n\n"; if(AllowFail) { return; } Error(); }
				CheckAlloc[pLeft] = true;
				m_Node[pLeft]->SetNode(NextPointer,pLeft,-1);
			} else	{
				assert(m_Node[pLeft]->m_viLink.size() > 2 && m_Node[pLeft]->m_viBranch.size() > 2);
//				cout << "\n\tDoing Left["<<m_Node[pLeft]->m_viLink[2]<<"]: " << NextPointer;
				m_Node[pLeft]->m_viLink[2] = NextPointer;
				m_Node[pLeft]->m_viBranch[2] = pLeft;
			}
			SetB(pLeft,TempBranch[0],true,true); m_viBranchLabels[pLeft] = TempLabel[0];
			if(pRight < m_iNoSeq)        {
				if(CheckAlloc[pRight] != false)     { cout << "\nNode " << pRight << " multiply allocated.\nUsually caused by same sequence occuring multiple times in tree\n\n"; if(AllowFail) { return; } Error(); }
				CheckAlloc[pRight] = true;
				m_Node[pRight]->SetNode(NextPointer,pRight, -1);
			} else        {
				assert(m_Node[pRight]->m_viLink.size() > 2 && m_Node[pRight]->m_viBranch.size() > 2);
//				cout << "\n\tDoing Right["<<m_Node[pRight]->m_viLink[2]<<"]: " << NextPointer;
				m_Node[pRight]->m_viLink[2] = NextPointer;
				m_Node[pRight]->m_viBranch[2] = pRight;  }
			SetB(pRight,TempBranch[1],true,true); m_viBranchLabels[pRight] = TempLabel[1];
			// Initialise the new Node
			if(CheckAlloc[NextPointer] != false)	{ cout << "\nNode " << NextPointer << " multiply allocated.\nUsually caused by same sequence occuring multiple times in tree\n\n"; if(AllowFail) { return; } Error(); }
			CheckAlloc[NextPointer] = true;
			// Create an internal node
//			cout << "\n\tAllocating internal node["<<NextPointer<<"]: L=" << pLeft << "; R=" << pRight << "; Par="<<Parent << flush;
			if(m_bRooted && i == countBra - 1) { m_Node[NextPointer]->SetNode(pLeft,pRight,pLeft,pRight,IntVal); }
			else { m_Node[NextPointer]->SetNode(pLeft,pRight,Parent,pLeft,pRight,NextPointer,IntVal); }
			NextPointer ++; // Move to next node
        }

		// For unrooted trees get the final 3 nodes of the tree
        if(!m_bRooted)	{
        	TempTree = TempTree.substr(1,(TempTree.find_last_of(")")-1));
        	Toks = Tokenise(TempTree,",");
        	if(Toks.size() != 3) { Error("When creating unrooted CTree: expected to see 3 nodes at base of tree, instead saw: " + TempTree + "\n\n"); }
        	FOR(i,3)	{
        		TempBranch[i] = SET_BRANCH;
        		GetBraDetails(Toks[i],&j,&TempBranch[i],&TempLabel[i]);
        		if(j < m_iNoSeq)	{	// If needed initialise leaf nodes
        			if(CheckAlloc[j] != false)	{ cout << "\nNode " << j << " multiply allocated.\nUsually caused by same sequence occuring multiple times in tree\n\n"; if(AllowFail) { return; } Error(); }
        			CheckAlloc[j] = true; m_Node[j]->SetNode(NextPointer,j, -1);
        		}
        		// Assign the branch lengths
        		SetB(j,TempBranch[i],true,true); m_viBranchLabels[j] = TempLabel[i];
        		switch(i) {
        			case 0:	pLeft = j; break;
        			case 1: pRight = j; break;
        			case 2: Parent = j; break;
        			default: cout << "\nError assigning branch lengths...\n\n"; if(AllowFail) { return; } Error();
        	}	}
        	if(m_Node[NextPointer] != NULL) { delete m_Node[NextPointer]; }
        	m_Node[NextPointer] = new CNode(pLeft,pRight,Parent,pLeft,pRight,Parent,IntVal); // Initialise the final node
        	if(m_Node[pLeft]->m_viLink.size() > 2)	{ m_Node[pLeft]->m_viLink[2] = NextPointer; }
        	if(m_Node[pRight]->m_viLink.size() > 2)	{ m_Node[pRight]->m_viLink[2] = NextPointer; }
        	if(m_Node[Parent]->m_viLink.size() > 2)	{ m_Node[Parent]->m_viLink[2] = NextPointer; }
//        	cout << "\nDone!" << flush;
        } else {
			// Do the rooted tree by adjusting the root node.
        	m_iRootNode = NextPointer -1;
        }
	}	else	{	// Otherwise for two species
		if(SubSeq == m_iNoSeq) {
			if(m_Node[0] != NULL) { delete m_Node[0]; }
			m_Node[0] = new CNode(1,0); m_Node[1] = new CNode(0,0);
			//cout << "\nTempTree: " << TempTree << " cf. " << TempBranch << flush;
			i=1; while(isdigit(TempTree[i])) { i++; }	// Do first side
			SetB(0,DoBranch(&TempTree,&i,&IntVal),true,true);
			while(TempTree[i]!=':' && TempTree[i] != ';')	{ i++; }
			if(TempTree[i] == ':') { AddB(0,DoBranch(&TempTree,&i,&IntVal),true,true); } // Get branch2
		}
		else {
			TempBranch[0] = TempBranch[1] = TempBranch[2] = 0.0;
			find_closest(&TempTree,&pLeft,&pRight,&Parent,NextPointer, TempBranch, TempLabel,&IntVal);
			if(m_Node[pLeft] != NULL) { delete m_Node[pLeft]; }
			m_Node[pLeft] = new CNode(pRight,0); m_Node[pRight] = new CNode(pLeft,0);
			SetB(0,TempBranch[0] + TempBranch[1],true,true);
			ReplaceBraLink(0,0,pRight); ReplaceBraLink(0,1,pLeft);
		}
	}
    FOR(i,NoNode()) { if(CheckAlloc[i] == false) { break; } } // Check if all nodes assigned
	if(i!=NoNode()) { CheckVar = false; FOR(i,m_iNoSeq) { if(NoLinks(i) > 0) {  m_iStartCalc = i; break; } } }
	BuildBraLinks(CheckVar);  // Get the list of branches links
	OrderNode(-1,true);
	m_bReady=true;	    // Tree is now ready!!!
	// Do some branch by branch things
	m_iNoOptBra = m_iNoBra;
	DEL_MEM(CheckAlloc);
	ValidateTree();
	m_bValid = true;
	// Remove branch labels if not used
	for(i=1;i<(int)m_viBranchLabels.size();i++) { if(m_viBranchLabels[i] != m_viBranchLabels[0]) { break; } }
	if(i == m_viBranchLabels.size()) { m_viBranchLabels.clear(); }
	// Sort branch labels so sensibly labelled
	temp_l.assign(m_viBranchLabels.begin(),m_viBranchLabels.end());
	temp_l.sort(); temp_l.unique();
	temp_l2.assign(temp_l.begin(),temp_l.end());
	FOR(i,(int)m_viBranchLabels.size()) { FOR(j,(int)temp_l2.size()) { if(temp_l2[j] == m_viBranchLabels[i]) { m_viBranchLabels[i] = j; break; } } }
	m_iNoBraLabels = (int)temp_l2.size();
}

// Create the branch links. Also works for sub trees
void CTree::BuildBraLinks(bool Verify)	{
	int i,j,k;
	vector <int> BranCount;
	// Initialise
	m_vBraLinks.clear();
	m_vBraLinks.assign(m_iNoBra,vector<int>(2,-1));
	BranCount.assign(m_iNoBra,0);
	// Get assignments
	FOR(i,m_iNoNode)	{
		FOR(j,(int)m_Node[i]->m_viBranch.size())	{
			k = m_Node[i]->m_viBranch[j];
			if(k == -1) { continue; }
			assert(k < m_iNoBra && k >= 0);
			if(m_vBraLinks[k][0] == -1)	{ m_vBraLinks[k][0] = i; BranCount[k]++;}
			else							{ m_vBraLinks[k][1] = i; BranCount[k]++;}
	}	}
	// Verify assignments
	if(Verify == true) {
		FOR(i,m_iNoBra) {
			if(BranCount[i] != 2) { cout << "\nError in BranCount["<<i<<"]: " << BranCount[i] << flush; cout << "\nTree: " << *this << flush; }
			assert(BranCount[i] == 2);
	} }
}

void CTree::BuildBranches(int NT, int NF)	{
	int i,Node;
	// If first call then sort it out
	if(NT == -1 && NF == -1)	{
		Node = (2*m_iNoSeq) - 3;
		// Go down the tree and do branches
		FOR(i,(int)m_Node[Node]->m_viLink.size()) { BuildBranches(m_Node[Node]->m_viLink[i],Node); }
		return;
	}
	// PostOrder Tree-Traversal
	FOR(i,(int)m_Node[NT]->m_viLink.size())	{
		if(m_Node[NT]->m_viLink[i] == NF || m_Node[NT]->m_viLink[i] == -1) { m_Node[NT]->m_viBranch[i] = NT; continue; }
		BuildBranches(m_Node[NT]->m_viLink[i],NT);
	}
	FOR(i,(int)m_Node[NF]->m_viLink.size())	{
		if(m_Node[NF]->m_viLink[i] == NT) { m_Node[NF]->m_viBranch[i] = NT; break; }
	}

}

// Find branch linking nodes i and j; returns -1 if doesn't exist
int CTree::FindBra(int N_i, int N_j)	{
    int i;
	FOR(i,m_iNoBra)	{
		if( (m_vBraLinks[i][0] == N_i && m_vBraLinks[i][1] == N_j) ||
			(m_vBraLinks[i][0] == N_j && m_vBraLinks[i][1] == N_i) )	{ return i; }
	}
	return -1;
}

// Create a set of labels for branches corresponding to their numbers
void CTree::CreateBranchLabels(bool Force) {
	int i;
	if(!Force) { FOR(i,(int)m_viBranchLabels.size()) { if(m_viBranchLabels[i] != m_viBranchLabels[0])  { Error("\nTrying to CTree::CreateBranchLabels when they already exist\n\n"); } } }
	if(m_viBranchLabels.empty()) { m_viBranchLabels.assign(m_iNoBra,0); }
	FOR(i,(int)m_viBranchLabels.size()) { m_viBranchLabels[i] = i; }
	m_iNoBraLabels = m_viBranchLabels.size();
}

// Destructor function
///////////////////////////////////////////////////
CTree::~CTree() {
#if DO_MEMORY_CHECK
	memory_check.CountCTree--;
#endif
	CleanTree();
}

// Allocate the memory
void CTree::GetMemory(int NoSeq)	{
	string Name;
	assert(m_bReady == false);
	// Adding an extra node to all trees because tree may be rooted
	int i;
	assert(m_Node.empty() && m_vBraLinks.empty() && m_vdBra.empty());
    if(m_bRooted) { m_iNoNode = (2*NoSeq)-1; m_iNoBra = (2*NoSeq)-2; m_iNoSeq = NoSeq; }
    else { m_iNoNode = (2*NoSeq)-2; m_iNoBra = (2*NoSeq)-3; m_iNoSeq = NoSeq; }
    m_vBraLinks.assign(m_iNoBra,vector<int>(2,-1));
    m_vdBra.assign(m_iNoBra,0.0);
	// Set up nodes
    m_Node.assign(m_iNoNode,NULL);
    for(i=0;i<(int)m_Node.size();i++) { m_Node[i] = NULL; m_Node[i] = new CNode; assert(m_Node[i] != NULL); }
}

// Clean tree function
void CTree::CleanTree() {
	int i;
    for(i=0;i<m_Node.size();i++) { if(m_Node[i] != NULL) { delete m_Node[i]; }}
    m_Node.clear();
    m_vBraLinks.clear();
	m_vdBra.clear();
	m_bReady = false; m_iOutBra = m_bOutName = false; m_bFastCalcOK = false; m_bOutLabel = false;
}

// Validate tree function
/////////////////////////////////////////////////
// 1. Checks all nodes have appropriate forwards and backwards relationships
// 2. Checks that all branches are linked where they think they are
bool CTree::ValidateTree(bool AllowExit)	{
#if TREE_DEBUG
	int Node,Link,i;
	FOR(Node,NoNode())	{
//		cout << "\nValidate: Node["<<Node<<"], link " << flush;
		FOR(Link,NoLinks(Node))	{
//			cout << "\t[" << Link << "] = " << NodeLink(Node,Link) << flush;
			if(NodeLink(Node,Link) != -1) {
				// Check the node looks back
				FOR(i,NoLinks(NodeLink(Node,Link)))	{ if(NodeLink(NodeLink(Node,Link),i) == Node) { break; } }
				if(i == NoLinks(NodeLink(Node,Link))) { OutDetail(cout,true); if(AllowExit) { Error("\nTree links invalid...\n\n"); } return false; }
			}
			if(NodeBra(Node,Link) != -1)	{
				// Check the branches and nodes match
				FOR(i,2) { if(m_vBraLinks[NodeBra(Node,Link)][i] == Node || m_vBraLinks[NodeBra(Node,Link)][i] == -1) { break; }
				if(i == 2) { OutDetail(cout,true); if(AllowExit) { Error("\nTree branches invalid...\n\n"); } return false; } }
	}	}	}
#endif
	return true;
}

//////////////////////////////////////////////////////////////////////
// Access functions
//////////////////////////////////////////////////////////////////////
vector <double> CTree::Branches()	{
	int i;
	vector <double> Br;
	FOR(i,NoBra()) { Br.push_back(B(i)); }
	return Br;
}
double CTree::TreeLength() {
	vector <double> Temp = Branches();
	return Sum(&Temp);
}
double CTree::B(int Branch)	{
	assert(InRange(Branch,0,(int)m_vdBra.size()));
	return m_vdBra[Branch];
}
double CTree::SetB(int Branch, double Value, bool Update, bool Rescale)	{
	assert(InRange(Branch,-1,(int)m_vdBra.size()));
	m_vdBra[Branch] = Value;
	return B(Branch);
}
double CTree::MulB(int Branch, double Value, bool Update, bool Rescale) {
	assert(InRange(Branch,0,(int)m_vdBra.size()));
	return SetB(Branch,Value * B(Branch), Update, Rescale); ;
}
double CTree::AddB(int Branch, double Value, bool Update, bool Rescale) {
	assert(InRange(Branch,0,(int)m_vdBra.size()));
	return SetB(Branch, Value + B(Branch), Update, Rescale);
}

double CTree::QuadB(int Branch)	{
	int i,j,NoBra = 1;
	double  TotalBra = B(Branch);	// Gets branch length (and checks the branch is valid
	FOR(i,2)	{
		FOR(j,(int)m_Node[m_vBraLinks[Branch][i]]->m_viLink.size())	{
			if(NodeBra(m_vBraLinks[Branch][i],j) == -1 || NodeLink(m_vBraLinks[Branch][i],j) == -1) { return -BIG_NUMBER; }
			if(NodeBra(m_vBraLinks[Branch][i],j) == Branch) { continue; }
			NoBra++;
			TotalBra += B(m_Node[m_vBraLinks[Branch][i]]->m_viBranch[j]);
	}	}
	return TotalBra / (double) NoBra;
}

void CTree::ReplaceNodeLink(int N,vector <int> L)	{
	assert(InRange(N,0,NoNode()));
	m_Node[N]->m_viLink  = L;
	m_bFastCalcOK = false;
}
void CTree::ReplaceNodeLinkElement(int N, int Element, int Val)	{
	assert(InRange(N,0,NoNode()));
	assert(InRange(Element,0,NoLinks(N)));
	m_Node[N]->m_viLink[Element] = Val;
	m_bFastCalcOK = false;
}
void CTree::ReplaceNodeBra(int N,vector <int> B)	{
	assert(InRange(N,0,NoNode()));
	m_Node[N]->m_viBranch = B;
	m_bFastCalcOK = false;
}
void CTree::ReplaceNodeBraElement(int N, int Element, int Val)	{
	assert(InRange(N,0,NoNode()));
	assert(InRange(Element,0,NoLinks(N)) && NoLinks(N) == m_Node[N]->m_viBranch.size());
	m_Node[N]->m_viBranch[Element] = Val;
	m_bFastCalcOK = false;
}

// Tree operator function
////////////////////////////////////////////////

void CTree::operator=(const CTree & Tree)	{
    int i;
	if(!Tree.m_bReady || &Tree == this) { return; }
	CleanTree();				// Clear memory
	m_bRooted = Tree.m_bRooted; m_iRootNode = Tree.m_iRootNode;
    GetMemory(Tree.m_iNoSeq);	// New memory
    // Transfer data
    for(i=0;i<m_iNoBra;i++)  {
		m_vdBra[i] = Tree.m_vdBra[i];
		m_vBraLinks[i][0] = Tree.m_vBraLinks[i][0];
		m_vBraLinks[i][1] = Tree.m_vBraLinks[i][1];
	}
	m_vsName = Tree.m_vsName;
	// Some values regarding models and optimisation
	m_bOutName = Tree.m_bOutName;
	m_bOutLabel = Tree.m_bOutLabel;
	m_iOutBra = Tree.m_iOutBra;
	m_iNoOptBra = Tree.m_iNoOptBra;
	m_viBranchLabels = Tree.m_viBranchLabels; m_iNoBraLabels = Tree.m_iNoBraLabels;
    for(i=0;i<m_iNoNode;i++) { *m_Node[i] = *Tree.m_Node[i]; }
	m_iStartCalc = Tree.m_iStartCalc;
	m_bReady = true; m_bOldTree = Tree.m_bOldTree; m_bFastCalcOK = false;
}

// Find closest function
///////////////////////////////////////////////////////////
// Find closest two nodes together in a string
// Arguements:
// 	tree = The string to get info from
// 	child1, child2 and parent are the linked nodes
// 	next_pointer = node for information to be allocated to
//////////////////////////////////////////////////////////////////
void CTree::find_closest(string * tree, int *c1, int *c2, int *p,int n_p, double bra[],int label[], int *IntVal) {
	int x,x2,y;
	string resolve, sub;
	x = tree->find_first_of("("); x2 = tree->find("(",x+1); y = tree->find_first_of(")") + 1;
	// Find the matching parentheses
	while(y > x2)	{
		x = x2; x2 = tree->find("(",x+1);
		if(x2 == string::npos) { x2 = (int)tree->size();}
	}
    // Get the stuff to resolve and replace it in the tree
	resolve = tree->substr(x,y-x);
	tree->replace(x,y-x,int_to_string(n_p +1));
	// Get the required values
	x2 = resolve.find_first_of(",");
	// 1. Do left hand side
	GetBraDetails(resolve.substr(resolve.find_first_of("(")+1,x2-1),c1,&bra[0],&label[0]);
	// 2. Do right hand side
	GetBraDetails(resolve.substr(x2+1),c2,&bra[1],&label[1]);
//	cout << "\nc1 = " << *c1 << ", c2 = " << *c2 << ", m_iNoNode: " << m_iNoNode;
//	cout << "\nNew tree: " << *tree;
    // Reset for next part of loop
    *p=n_p;
	assert(InRange(*c1,0,m_iNoNode) && InRange(*c2,0,m_iNoNode));
}

// Function for getting branch details from a string of form Number:length#label
// Assumes bra and node are correctly initialised
void CTree::GetBraDetails(string s, int *node, double *bra, int *label)	{
	int c,d,i,j;
	string sub;
	*label = 0; *bra = SET_BRANCH;
	// Organise '(' and ')' in strings; strictly one of each allowed
	i = s.find("("); j = s.find(")");
	if(i!= string::npos) { s.replace(i,1,""); } i = s.find("("); if(j!= string::npos) { s.replace(j,1,""); }
	i = s.find("("); j = s.find(")");
	if(i != string::npos || j != string::npos) { Error("\nCTree::GetBraDetails(...)\nTrying to resolve string with more than one parenthesis: " + s + " \n"); }
//	cout << "\nGetBraDetails: " << s;
	i = s.find(":"); if(i==string::npos) { i = s.size(); }
	j = s.find(":"); if(j==string::npos) { j = s.size(); }
	*node = atoi(s.substr(0,min(i,j)).c_str()) - TREE_START;
	for(c=min(i,j);c<(int)s.size();c++) {
		if(s[c] == ':') {
			for(d=c+1;d<(int)s.size();d++)  {
				if(!(isdigit(s[d]) || s[d] == '.' || toupper(s[d]) == 'E' || s[d] == '-')) { break; } }
			*bra = atof(s.substr(c+1,d-c-1).c_str());
		} else if (s[c] == '#') {
			for(d=c+1;d<(int)s.size();d++)  { if(!(isdigit(s[d]) || s[d] == '.')) { break; } }
			*label = atoi(s.substr(c+1,d-c-1).c_str());
		}
	}
//	cout << "\nNode: " << *node << " bra: " << *bra << " label: " << *label;
}

double CTree::DoBranch(string *resolve,int *pos, int *IntVal)	{
	int i = *pos,j;
	double ret_val = SET_BRANCH;
	// If there is obviously no branch return
	if(i+1 >= (int)resolve->size()) { return ret_val; }
	if(resolve->at(i) ==':') {
		i++;
		ret_val=atof(resolve->c_str()+i);
		// Skips the digits involved for this
		while((resolve->at(i)=='-')||(isdigit(resolve->at(i))) ||(resolve->at(i)=='.')||(resolve->at(i)=='e')) { i++;}
	}
	// Get the IntVal if it exists
	// Requirement is that it falls between ')' and branch length's ':'
	if(resolve->at(i) ==')' && resolve->at(i+1) != ':' && i+1 < (int)resolve->size()) {
		j = i; i++;
		while(resolve->at(i) != ':' && resolve->at(i) != ',' && resolve->at(i) != ')' && resolve->at(i) != '(' && resolve->at(i) != ';') { i++; }
		if(resolve->at(i)==':' || resolve->at(i) == ';') { *IntVal = atoi(resolve->c_str()+j+1); }
		i = j;
	}
	*pos = i;
	if(ret_val < 0) { ret_val = 0; }
	return ret_val;
}

// Removes a node from a tree that has already had one of its links removed
int CTree::RemoveBraNode(int N)	{
	int i,j, count;
	double BVal;
	vector <int> Link, Bra;
	// Initialise and check entry conditions
	assert(N >= NoSeq());
	count = 0;
	FOR(i,NoLinks(N)) {
		if(NodeLink(N,i) == -1) { count++; continue; }
		Link.push_back(NodeLink(N,i)); Bra.push_back(NodeBra(N,i));
		ReplaceNodeLinkElement(N,i,-1); ReplaceNodeBraElement(N,i,-1);
	}
	assert(count == 1 && Link.size() == 2 && Bra.size() == 2);
	// Do the branch lengths
	BVal = B(Bra[0])+B(Bra[1]);
	FOR(i,2) { SetB(Bra[i],BVal,true,true); }
	// Now remove the nodes
	FOR(i,2) {
		FOR(j,NoLinks(Link[i])) { if(NodeLink(Link[i],j) == N) { break; } }
		assert(j != NoLinks(Link[i]));
		ReplaceNodeLinkElement(Link[i],j,Link[FlipBin(i)]);
		ReplaceNodeBraElement(Link[i],j,Bra[0]);
		m_vBraLinks[Bra[0]][i] = Link[i];
	}
	OrderNode(Link[0]); OrderNode(Link[1]);
	m_vBraLinks[Bra[1]][0] = m_vBraLinks[Bra[1]][1] = -1;
	return Bra[0];
}

// Cuts a branch
//  returns the branch on priority subtree where they used to be attached -- i.e. the blank node
// -- Subtrees are maintained; Link specifies which one is given priority
int CTree::CutBranch(int B, int Link)	{
	int i,j,N,F, RetVal = -1;
	// Remove facing nodes
	FOR(i,2) {
		N = BraLink(B,i);	// Get node to
		F = BraLink(B,FlipBin(i));
		FOR(j,NoLinks(N)) { if(NodeLink(N,j) == F) { break; } }
		assert(j!=NoLinks(N));
		m_Node[N]->m_viLink[j] = -1; m_Node[N]->m_viBranch[j] = -1;
		if(i != Link) { RetVal = RemoveBraNode(N); }
		else { OrderNode(N); }
	}
	// Clean up the branch
	m_vBraLinks[B][0] = m_vBraLinks[B][1] = -1;
	return RetVal;
}

// Cut Sequence function: removes a sequence from the tree
///////////////////////////////////////////////////////////////
// Needs to store:
//	2 x branch numbers (leading to leaf node; trivial + spare branch)
//	1 x node number (node removed from tree)
// This will be held in the node (which is ready for addition to tree)
//
// Return value:
//	Spare node number
//
// Assumptions:
// The branch leading to a leaf node has the # of the leaf node

int CTree::RemoveLeafNode(int RemNode)	{
	/* Removes a Node from the tree
			returns 0 if all okay, else 1	*/
	int i,j,k,m,links[2],base_node;double total;
	// Error checking for value passed
	if(m_iNoSeq < 3) {
		cout << "\nAttempted to remove a sequence when only two are left";
		return 1;
	}
	if(RemNode >= m_iNoSeq)		{
		cout << "\nAttempted to remove an internal node";
		return 1;
	}
	// Section to remove the corresponding section from the tree
	// Deal with when removing a sequence from a three species tree
	if(m_iNoSeq == 3)	{
		if(RemNode == 0)		{ m_Node[0] = m_Node[2]; total = B(1) + B(2); }
		else if(RemNode == 1)	{ m_Node[1] = m_Node[2]; total = B(0) + B(2); }
		else	{ total = B(0) + B(1); }
		m_Node[0]->m_viLink[0] = 1; m_Node[0]->m_viBranch[0] = 0;
		m_Node[1]->m_viLink[1] = 0; m_Node[0]->m_viBranch[0] = 0;
		m_vBraLinks[0][0] = 0; m_vBraLinks[0][1] = 1;
		SetB(0,total,true,true);
		m_iNoSeq = 2;
		m_bFastCalcOK = false;
		return 0;
	}
	// Remove a sequence from the tree
	// Get base_node and links
	base_node = m_Node[RemNode]->m_viLink[0];
	assert(m_Node[base_node]->m_NodeType == branch);	j=0;
	for(total=0.0,i=0;i<3;i++)	{
		if(m_Node[base_node]->m_viLink[i] == RemNode) { continue; }
		links[j++] = m_Node[base_node]->m_viLink[i];
		if(j==2) { // Set the disconnected branch link to null
			FOR(m,2) {
				if(m_vBraLinks[ m_Node[base_node]->m_viBranch[i] ][m] != links[1]) { m_vBraLinks[ m_Node[base_node]->m_viBranch[i] ][m] = -1; }
			}	}
		total += B(m_Node[base_node]->m_viBranch[i]);
	}
	if(links[0] > links[1]) { i = links[0]; links[0] = links[1]; links[1] = i; }
	// Reassign the correct links and so on
	// Important do not change iOption(sequence node removed) or base_node(node removed)
	k = 0;  // the altered nodes will point to
	for(k=-1,i=0;i<2;i++)	{
		if(m_Node[links[i]]->m_NodeType == branch)	{
			for(j=0;j<3;j++) {
				if(m_Node[links[i]]->m_viLink[j] == base_node)	{	// Do rearrangements
					if(i == 0)	{
						m_Node[links[0]]->m_viLink[j] = links[1];
						k = m_Node[links[0]]->m_viBranch[j];
						SetB(k,total,true,true);
						FOR(m,2)	{ // Reassign branch links
							if(m_vBraLinks[k][m] == base_node) { m_vBraLinks[k][m] = links[1]; break; }
					}	}
					else		{
						assert(k>=0);
						m_Node[links[1]]->m_viLink[j] = links[0];
						m_Node[links[1]]->m_viBranch[j] = k;
		}	}	}	}
		else	{
			if(i==0) {
				m_Node[links[0]]->m_viLink[0] = links[1];
				k = m_Node[links[0]]->m_viBranch[0];
				SetB(k,total,true,true);
				FOR(m,2)	{ // Reassign branch links
					if(m_vBraLinks[k][m] == base_node) { m_vBraLinks[k][m] = links[1]; break; }
			}	}
			else	{
				m_Node[links[1]]->m_viLink[0] = links[0];
				// This shouldn't have to be here...
				SetB(m_Node[links[1]]->m_viBranch[0],total,true,true);
				m_Node[links[1]]->m_viBranch[0] = k;
	}	}	}
	// Tidy up the basenode
	FOR(i,3)	{
		if(m_Node[base_node]->m_viLink[i] == RemNode) { continue; }
		if(m_Node[base_node]->m_viLink[i] == links[0]) { // This is the branch thats removed
			m_Node[base_node]->m_viBranch[i] = -1;
		}
		m_Node[base_node]->m_viLink[i] = -1;
	}
	OrderNode();
	BuildBraLinks(false);
	m_bFastCalcOK = false;
	ValidateTree();
	return 0;
}

// AddSeq
// ------
// For a sequence to be added it requires a bit of tree structure associated with it:
// 1) two nodes: one associated with the sequence and one internal
// 2) two branches: linking leaf and internal node, and a spare.
// BranchProp = relative position in branch (from low to high node nums)
// This must be prepared before passing to the function
//
// Returns -1 if the branch is the same as one of the ones in the cut section
int CTree::AddSeq(int SNum,int Bra, double Prop)	{
	int i,j,braN, braB, leafB,inN1, inN2,okay1, okay2;;
	double Len = (1.0 - Prop) * B(Bra);
	// Check entry conditions w.r.t. tree and initialise some variables
	assert(SNum < m_iNoSeq && m_Node[SNum] != NULL && SNum < m_iNoNode);
	assert(m_Node[SNum]->m_NodeType == leaf && (!m_Node[SNum]->m_viLink.empty()) && (!m_Node[SNum]->m_viBranch.empty()));
	braN = m_Node[SNum]->m_viLink[0]; leafB = m_Node[SNum]->m_viBranch[0];
	assert(braN < m_iNoNode && m_Node[braN]->m_NodeType == branch && m_Node[braN]->m_viLink.size() == 3 && m_Node[braN]->m_viBranch.size() == 3);
	braB = m_Node[braN]->m_viBranch[1];
	assert(m_Node[braN]->m_viLink[0] == SNum && m_Node[braN]->m_viLink[1] == -1 && m_Node[braN]->m_viLink[2] == -1);
	assert(m_Node[braN]->m_viBranch[0] == leafB);
	assert(m_Node[braN]->m_viBranch[0] == leafB && m_Node[braN]->m_viBranch[1] != -1 && m_Node[braN]->m_viBranch[2] == -1);
	inN1 = m_vBraLinks[Bra][0]; inN2 = m_vBraLinks[Bra][1];
	assert(IsNode(inN1) && IsNode(inN2));
	if(inN1 > inN2)  { i = inN1; m_vBraLinks[Bra][0] = inN1 = inN2; m_vBraLinks[Bra][1] = inN2 = i; }

	// Adjust some branch lengths
	MulB(Bra,Prop,true,true);
	SetB(braB,Len,true,true);
	// Adjust the BraN
	m_Node[braN]->m_viLink[1] = inN2; m_Node[braN]->m_viLink[2] = inN1;
	m_Node[braN]->m_viBranch[2] = Bra;
	// Adjust inN1
	FOR(i,(int)m_Node[inN1]->m_viLink.size()) { if(m_Node[inN1]->m_viLink[i] == inN2) { break; } } assert(i!=m_Node[inN1]->m_viLink.size());
	m_Node[inN1]->m_viLink[i] = braN;
	FOR(j,2) { if(m_vBraLinks[Bra][j] != inN1) { break; } } assert(j!=2);
	m_vBraLinks[Bra][j] = braN;
	// Adjust inN2
	FOR(i,(int)m_Node[inN2]->m_viLink.size()) { if(m_Node[inN2]->m_viLink[i] == inN1) { break; } } assert(i!=m_Node[inN2]->m_viLink.size());
	m_Node[inN2]->m_viLink[i] = braN;
	m_Node[inN2]->m_viBranch[i] = braB;
	FOR(j,2) { if(m_vBraLinks[braB][j] != braN) { break; } } assert(j!=2);
	m_vBraLinks[braB][j] = inN2;
	// Check node have been assigned correctly
	okay1 = okay2 = 0;
	assert(m_vBraLinks[leafB][0] == SNum && m_vBraLinks[leafB][1] == braN);
	assert(m_vBraLinks[Bra][0] == inN1 && m_vBraLinks[Bra][1] == braN);
	okay1 = 0; FOR(i,2) { if(m_vBraLinks[braB][i] == inN2) { okay1++; } else if(m_vBraLinks[braB][i] == braN) { okay1++; } }
	assert(okay1 == 2);
	okay1 = okay2 = 0;
	FOR(i,3) { 		assert(IsNode(m_Node[braN]->m_viLink[i]) && IsBra(m_Node[braN]->m_viBranch[i])); }
	FOR(i,(int)m_Node[inN1]->m_viLink.size())	{
		if(m_Node[inN1]->m_viLink[i] == braN) { okay1++; }
		assert(IsNode(m_Node[inN1]->m_viLink[i]) && IsBra(m_Node[inN1]->m_viBranch[i]));
	}
	FOR(i,(int)m_Node[inN2]->m_viLink.size())	{
		if(m_Node[inN2]->m_viLink[i] == braN) { okay2++; }
		assert(IsNode(m_Node[inN2]->m_viBranch[i]) && IsBra(m_Node[inN2]->m_viBranch[i]));

	}
	assert(okay1 = 1 && okay2 == 1);
	OrderNode();
	m_bFastCalcOK = false;
	ValidateTree();
	return 1;
}

bool CTree::IsNode(int Node) {
	if(Node >= 0 && Node < m_iNoNode) { return true; }
	return false;
}

bool CTree::IsBra(int Bra)	{
	if(Bra >= 0 && Bra < m_iNoBra) { return true; }
	return false;
}

bool CTree::GoodBra(int Bra)	{
	int i;
	FOR(i,2)	{
		if(m_vBraLinks[Bra][i] == -1) { return false; }
		if(GoodNode(m_vBraLinks[Bra][i]) == false) { return false; }
	}
	return true;
}

bool CTree::GoodNode(int Node)	{
	int NodNum = Node;
	if(Node < m_iNoSeq) { NodNum = m_Node[Node]->m_viLink[0]; }
	if(m_Node[NodNum]->m_viLink[2] == -1) { return false; }
	return true;
}

int CTree::NoLeafLink(int N)	{
	int i,count = 0;
	FOR(i,NoLinks(N)) { if(NodeLink(N,i) < m_iNoSeq) { count++; } }
	return count;
}

bool CTree::IsCutTree()	{
	int i;
	FOR(i,m_iNoBra)	{
		if(m_vBraLinks[i][0] == -1 || m_vBraLinks[i][1] == -1) { return true; }
	}
	return false;
}

int CTree::GetStart(bool replace)	{
	int i;
	FOR(i,m_iNoSeq) {
		if(m_Node[i] == NULL) { continue; }
		if((int)m_Node[i]->m_viLink.size() !=1 ) { continue; }
		if(!InRange(m_Node[i]->m_viLink[0],0,m_iNoNode)) { continue; }
		if((int)m_Node[m_Node[i]->m_viLink[0]]->m_viLink.size() != 3) { continue; }
		if(m_Node[m_Node[i]->m_viLink[0]]->m_viLink[2]!=-1) { break; }
	}
	assert(i != m_iNoSeq);
	if(replace == true) { m_iStartCalc = i; }
	return i;
}

double CTree::GetTreeLength(bool first, int NTo, int NFr)	{
	int i,NodeNum;
	static double length = 0;
	// Initialise function
	if(first == true) {
		length = 0;
		if(m_iNoSeq == 2) { return B(0); }	// Always just the first branch for 2 sequences
		return GetTreeLength(false,GetStart(false),-1);
	}
	assert(NTo >= 0 && NTo < m_iNoNode);
	if(m_Node[NTo]->m_NodeType == branch) { NodeNum = 3; } else { NodeNum = 1; }
	FOR(i,NodeNum)	{
		if(m_Node[NTo]->m_viLink[i] == NFr) { continue; }
		length += B(m_Node[NTo]->m_viBranch[i]);
		GetTreeLength(false,m_Node[NTo]->m_viLink[i],NTo);
	}
	return length;
}

int CTree::BranchSets(int BranchNum, vector <int> &Left, vector <int> &Right)	{
	vector <int> temp, temp2;
	assert(BranchNum < m_iNoBra);
	BuildBraLinks(false);
	GetBraSets(m_vBraLinks[BranchNum][0], m_vBraLinks[BranchNum][1],Left, true);
	GetBraSets(m_vBraLinks[BranchNum][1], m_vBraLinks[BranchNum][0],Right, true);
	// Sort the nodes
	if((int)Left.size() > 1) { sort(Left.begin(),Left.end()); }   else if(Left.empty()) { return -1; } else if(Left[0] == -1) { return -1; }
	if((int)Right.size() > 1) { sort(Right.begin(),Right.end()); } else if(Right.empty()) { return -1; } else if(Right[0] == -1) { return -1; }
	// Ensure the one with sequence 0 in is in the left
	if(Left[0] != 0) { temp = Left; Left = Right; Right = temp; }
	if(Left.size() > Right.size()) { return (int)Left.size(); }
	temp.clear(); temp2.clear();
	return (int) Right.size();
}

void CTree::GetBraSets(int NTo, int NFr, vector <int> &List, bool First)	{
	int i;
	if(NTo == -1) { return; }	// Avoid empty links in the tree
	if(First == true) { List.clear(); if(NTo < m_iNoSeq) { List.push_back(NTo); } }
	FOR(i,(int)m_Node[NTo]->m_viLink.size())	{
		// If already seen then continue
		if(m_Node[NTo]->m_viLink[i] == NFr || m_Node[NTo]->m_viLink[i] == -1) { continue; }
		// If an internal node then descend
		else if(m_Node[NTo]->m_viLink[i] >= m_iNoSeq)	{ GetBraSets(m_Node[NTo]->m_viLink[i],NTo,List,false); }
		// If an external node then store
		else { List.push_back(m_Node[NTo]->m_viLink[i]); }
}	}

// Order the links in a node
void CTree::OrderNode(int NodeNum, bool DoBraToo)	{
	int i,j,k,tempL,tempB, Valj, Valk;
	for(i=m_iNoSeq; i<m_iNoNode; i++)	{
		if(NodeNum == i || NodeNum == -1)	{
			// Sort nodes links and branches
			assert(m_Node[i]->m_viLink.size() == m_Node[i]->m_viBranch.size());
			FOR(j,(int)m_Node[i]->m_viLink.size())	{
				for(k=j;k<(int)m_Node[i]->m_viLink.size();k++)	{
					Valj = m_Node[i]->m_viLink[j]; Valk = m_Node[i]->m_viLink[k];
					if(Valj == -1)  { Valj = BIG_NUMBER; }
					if(Valk == -1)  { Valk = BIG_NUMBER; }
					if(Valj> Valk) {
						tempL = m_Node[i]->m_viLink[j]; tempB = m_Node[i]->m_viBranch[j];
						m_Node[i]->m_viLink[j] = m_Node[i]->m_viLink[k]; m_Node[i]->m_viBranch[j] = m_Node[i]->m_viBranch[k];
						m_Node[i]->m_viLink[k] = tempL; m_Node[i]->m_viBranch[k] = tempB;
			}	}	}
			if(m_Node[i]->m_viLink.size() == 3)	{
				if(m_Node[i]->m_viLink[1] == -1 && m_Node[i]->m_viLink[2] == -1)	{
					if(m_Node[i]->m_viBranch[2] > m_Node[i]->m_viBranch[1]) {
						tempB = m_Node[i]->m_viBranch[1];
						m_Node[i]->m_viBranch[1] = m_Node[i]->m_viBranch[2];
						m_Node[i]->m_viBranch[2] = tempB;
	}	}	}	}	}
	if(DoBraToo == true) {	// If need to do branches
		FOR(i,m_iNoBra)	{
			if( ( m_vBraLinks[i][0] > m_vBraLinks[i][1] && m_vBraLinks[i][1] != -1) || m_vBraLinks[i][0] == -1)	{
				j = m_vBraLinks[i][0];
				m_vBraLinks[i][0] = m_vBraLinks[i][1];
				m_vBraLinks[i][1] = j;
}	}	}	}

///////////////////////////////////////////////////////////////////
// Tree rooting and unrooting options
void CTree::Unroot()	{
	if(!IsRooted()) { return; }
	if(!InRange(m_iRootNode,0,m_iNoNode)) { Error("\nTrying to CTree::Unroot when tree not ready\n\n"); }
	assert(m_Node[m_iRootNode]->m_NodeType == root);
	int i, branch = NodeBra(m_iRootNode,0), pLeft = NodeLink(m_iRootNode,0), pRight = NodeLink(m_iRootNode,1);
	int braLeft = NodeBra(m_iRootNode,0), braRight = NodeBra(m_iRootNode,1);
	// Correct the branch length and remove the right hand branch
	m_vdBra[pLeft] = m_vdBra[braLeft] + m_vdBra[braRight];
	m_vdBra.erase(m_vdBra.begin() + pRight);
	// Fix left node
	FOR(i,NoLinks(pLeft)) {
		if(NodeLink(pLeft,i) == m_iRootNode) {
			ReplaceNodeLinkElement(pLeft,i,pRight);
			ReplaceNodeBraElement(pLeft,i,branch);
			break;
	}	}
	assert(i != NoLinks(pLeft));
	// Fix right node
	FOR(i,NoLinks(pRight)) {
		if(NodeLink(pRight,i) == m_iRootNode) {
			ReplaceNodeLinkElement(pRight,i,pLeft);
			ReplaceNodeBraElement(pRight,i,branch);
			break;
	}	}
	assert(i != NoLinks(pRight));
	// Now remove the redundant node
	delete m_Node[m_iRootNode];
	m_Node.erase(m_Node.begin() + m_iRootNode);
	// Finish up the other stuff
	m_iRootNode = -1; m_bRooted = false; m_iNoNode--; m_iNoBra--;
	BuildBraLinks(true); // Get the list of branches links
}

///////////////////////////////////////////////////////////////////
// Function that adds a midpoint root
void CTree::MidpointRoot() {
	double max_val = 0;
	int max_count = -1, count = 0;
	vector <double> pw = GetTreePW();
	for(auto & d : pw) { if(d >= max_val) { max_val = d; max_count = count;  } count++; }
	assert(max_count != -1);
	vector <int> no = GetNodePath(max_count % NoSeq(),max_count / NoSeq());
	max_val /= 2.0;
	for(int i = 0 ; i < no.size() - 1; i++) {
		int br = FindBra(no[i],no[i+1]);
		if(max_val - B(br) <= 0) {
			AddRoot(no[i],max_val,no[i+1]);
			break;
		}
		max_val -= B(br);
	}
}

void CTree::AddRoot(int NodeLeft, double BrLeft, int NodeRight) {
	int BrNumLeft = FindBra(NodeLeft,NodeRight), BrNumRight = m_iNoBra;
	// Check some conditions
	assert(BrLeft < B(BrNumLeft)); assert(!m_bRooted); assert(m_iRootNode == -1);
	// Initialise the other stuff
	m_iRootNode = m_Node.size(); m_bRooted = true; m_iNoNode++; m_iNoBra++;
	// Now create the nodes and branches
	CNode *NewNode; NewNode = new CNode;
	NewNode->SetNode(NodeLeft,NodeRight,BrNumLeft,BrNumRight,-1);
	m_Node.push_back(NewNode);
	m_vdBra.push_back(B(BrNumLeft) - BrLeft);
	m_vdBra[BrNumLeft] = BrLeft;
	// Change the attachments in the existing nodes
	for(auto &val : m_Node[NodeLeft]->m_viLink) { if(val == NodeRight) { val = m_iRootNode; } }
	for(auto &val : m_Node[NodeRight]->m_viLink) { if(val == NodeLeft) { val = m_iRootNode; } }
	for(auto &val : m_Node[NodeRight]->m_viBranch) { if(val == BrNumLeft) { val = BrNumRight; } }
	BuildBraLinks(true); // Get the list of branches links
}

vector <int> CTree::GetNodePath(int To, int From, vector <int> current, bool First) {
	static bool found;
	static int x, y;
	// Initialisation
	if(First) {
		x = To; y = From; From = -1;
		assert(InRange(x,0,NoNode()) && InRange(y,0,NoNode()));
		found = false;
	}
	if(To == y) { current.push_back(To); found = true; }
	// And traverse
	for(auto link : m_Node[To]->m_viLink) {
		if(link == From || found) { continue; }
		current = GetNodePath(link,To,current,false);
	}
	if(found && From != -1) { current.push_back(From); }
	return current;
}


vector <int> CTree::GetBranchPath(int To, int From, vector <int> current, bool First) {
	vector <int> nodes = GetNodePath(To,From);
	current.clear();
	for(int i = 0 ; i < nodes.size() - 1; i++) {
		current.push_back(FindBra(nodes[i],nodes[i+1]));
	}
	return current;
}

///////////////////////////////////////////////////////////////////
// Functions associated with splits on a tree
vector <SSplit> CTree::BuildSplits()	{
	int i;
	int rootCount = 0;
	m_vSplits.clear();
	FOR(i,NoBra()) {
		// Skip cases where the root leads to a trivial split
		if((BraLink(i,0) == Root() && BraLink(i,1) < NoSeq()) || (BraLink(i,1) == Root() && BraLink(i,0) < NoSeq())) {
			continue;
		}
		m_vSplits.push_back(CalculateSplit(i));
		if(BraLink(i,0) == Root() || BraLink(i,1) == Root()) {
			if(rootCount ++ == 0) {
				m_vSplits[i].rootLeft = true;
			} else {
				m_vSplits[i].rootRight = true;
			}
		}

	}
	assert(!(IsRooted() && rootCount == 0));
	return m_vSplits;
}

SSplit CTree::CalculateSplit(int Bra) {
	SSplit RetSplit;
	// Check entry conditions
	assert(InRange(Bra,0,m_iNoBra));
	// Get the splits
	RetSplit.BrLabel = Bra;
	BranchSets(Bra,RetSplit.Left,RetSplit.Right);
	if(!m_bRooted) {
		return RetSplit;
	}
//	if (Bra == m_Node[m_iRootNode]->m_viBranch[0] || Bra == m_Node[m_iRootNode]->m_viBranch[1]) { return RetSplit; } // Nothing for root
	SSplit root;
	BranchSets(m_Node[m_iRootNode]->m_viBranch[0], root.Left, root.Right);
	bool inLeft = false, inRight = false;
	// Test whether the root is attached to the left of RetSplit
	for(auto left : root.Left) {
		if(find(RetSplit.Left.begin(),RetSplit.Left.end(),left) != RetSplit.Left.end()) {inLeft = true; break;}
	}
	for(auto right : root.Right) {
		if(find(RetSplit.Left.begin(),RetSplit.Left.end(),right) != RetSplit.Left.end()) {inRight = true; break;}
	}
	if (inLeft && inRight) { RetSplit.rootLeft = true; } else { RetSplit.rootRight = true; }
	assert(!(RetSplit.rootLeft && RetSplit.rootRight));
	return RetSplit;
}

SSplit CTree::GetSplit(int Bra, bool forceRebuild) {
	if(m_vSplits.empty() || forceRebuild) { BuildSplits(); }
	if(InRange(Bra,0,(int)m_vSplits.size())) { return m_vSplits[Bra]; }
	return SSplit();			// Returns blank when there's a rooted tree and the root is on a terminal branch. Don't know if I like this behaviour
}

void CTree::OutSplits(ostream &os) {
	int i;
	BuildSplits();
	FOR(i,m_vSplits.size()) {
		cout << "\n\t" << m_vSplits[i].BrLabel << ":\t" << m_vSplits[i].Left << "\t" << m_vSplits[i].Right;
	}

}

///////////////////////////////////////////////////////////////////
// Function to calculate Robinson-Foulds distance between trees
int CTree::GetRFDist(CTree &Tree)	{
	int i,j, Dist;
	vector <vector <int> > L1, L2;
	vector <int> L,R;
	// Check some entry conditions
	if(NoSeq() != Tree.NoSeq()) { return -1; }
	// Get this trees sets
	FOR(i,NoBra()) { BranchSets(i,L,R); L1.push_back(L); Tree.BranchSets(i,L,R); L2.push_back(L); }
	// Get the distance
	Dist = NoBra();
	FOR(i,NoBra()) {
		FOR(j,NoBra()) { if(Compare(&L1[i],&L2[j]) == true) { Dist--; break; } }
	}
	return Dist;
}

////////////////////////////////////////////////////////////////////////////
// Function to check whether a subtree is compatible with the current tree
bool CTree::IsCompatible(CTree &SubTree) {
	// Check input information
	assert(SubTree.NoSeq() <= NoSeq());		// Check numbers of sequences are compatible
	// Build the sets of splits
	BuildSplits(); SubTree.BuildSplits();
	return SplitsCompatible(m_vSplits,NoSeq(),SubTree.m_vSplits,SubTree.NoSeq());
}

bool SplitsCompatible(vector <SSplit> Split1, int S1_seq, vector <SSplit> Split2, int S2_seq) {
	int i,j;
	// Check some entry conditions
	assert((int)Split1.size() >= (int)Split2.size());
	assert((int)Split1.size() == (S1_seq * 2) -3);
	assert((int)Split2.size() == (S2_seq * 2) -3);
//	cout << "\n>>>>>>>> INTO SplitsCompatible <<<<<<<<<<<";
	// Compare the subtree's splits to the full tree one-by-one to check they're compatible
	FOR(i,(int)Split2.size()) {
		if(Split2[i].Left.empty() || (int)Split2[i].Left.size() < 2 || (int)Split2[i].Right.size() < 2 ) { continue; }	// Skip empty splits or trivial splits
//		cout << "\nTesting SubTree["<<i<<"]: " << Split2[i].Left << " | " << Split2[i].Right;
		// Splits are compatible providing there one split from *this matches one split from SubTree
		// If not, then trees are incompatible and return false
		FOR(j,Split1.size()) {
			if((int)Split1[j].Left.size() < 2 || (int)Split1[j].Right.size() < 2) { continue; }
//			cout << "\n\tcf ["<<j<<"]: " << Split1[j].Left << " | " << Split1[j].Right;
			if(CompareSplit(Split1[j],Split2[i])) { break; }
		}
		// If loop runs to completion then no match is found
		if(j==(int)Split1.size()) { return false; }
	}
	// If all match then return true
	return true;
}

/////////////////////////// Compare Split /////////////////////
// Functions where S2 can have a subset of sequences from S1
// Assumes sequences in splits are in correct order (low to high); THIS IS NOT CHECKED!
bool CompareSplit(SSplit S1, SSplit S2) {
	int i, count = 0;
	bool MatchLeft = false;
	// Weak error checking
	assert((int)S1.Left.size() + (int)S1.Right.size() >= (int)S2.Left.size() + (int)S2.Right.size());
	// Compares S2.Left to S1.Left and S1.Right
	// 1. S2.Left vs S1.Left
	count = 0;
	FOR(i,(int)S1.Left.size()) {
		if(S1.Left[i] > S2.Left[count]) { break; }		// Always in order
		if(S1.Left[i] == S2.Left[count]) { count++; if(count == (int)S2.Left.size()) { MatchLeft = true; break; } }	// Match
	}
	// 2. S2.Left to S1. Right
	if(!MatchLeft) {
		count = 0;
		FOR(i,(int)S1.Right.size()) {
			if(S1.Right[i] > S2.Left[count]) { break; } 	// Always in order
			if(S1.Right[i] == S2.Left[count]) { count++; if(count == (int)S2.Left.size()) { break; } }	// Match
		}
		if(count != (int) S2.Left.size()) { return false; } // Neither Left nor right match
	}
	// Now checks right
	if(MatchLeft) {
		count = 0;
		FOR(i,(int)S1.Right.size())	{
			if(S1.Right[i] > S2.Right[count]) 	{ break; } 		// Always in order
			if(S1.Right[i] == S2.Right[count])	{ count++; if(count == (int)S2.Right.size()) { return true; } } 	// Match
		}
	} else {
		// 2. S2.Right vs S1.Left
		count = 0;
		FOR(i,(int)S1.Left.size()) {
			if(S1.Left[i] > S2.Right[count])	{ return false; } 	// Always in order
			if(S1.Left[i] == S2.Right[count])	{ count++; if(count == (int)S2.Right.size()) { return true; } }		// Match
		}
	}
	return false;
}

// ofstream operator for tree
/////////////////////////////////////////////////////
// NOTE: The code for this section is hideous
// 	 and badly described.  Can't be bothered
// 	 to debug however 'cause speed not issue
// 	 and it works
/////////////////////////////////////////////////////

ostream& operator<<(ostream& os, CTree &Tree)		{
	Tree.ValidateTree();
	if(Tree.m_iNoSeq == 2) {
		os << "(";
		if(Tree.m_bOutName && !Tree.m_vsName.empty()) { os << Tree.m_vsName[0]; } else { os << "1"; }
		if(Tree.m_iOutBra == 1) { os << ":" << Tree.B(0); }
		os << ",";
		if(Tree.m_bOutName && !Tree.m_vsName.empty()) { os << Tree.m_vsName[1]; } else { os << "2"; }
		if(Tree.m_iOutBra == 1) { os << ":0.0"; }
		os << ")";
		return os;
	}
	os << "(";
	if(Tree.IsCutTree())	{
		if(Tree.NoLinks(Tree.StartCalc()) == 0) { Error("Trying to start calc from bad node\n"); }
		Tree.OutNode(-1,Tree.m_Node[Tree.StartCalc()]->m_viLink[0],os);
	} else					{
		if(Tree.m_bRooted) { Tree.OutNode(-1,Tree.m_iRootNode,os); }
		else { Tree.OutNode(-1,Tree.m_Node[0]->m_viLink[0],os); }
		}
	os << ");";
	return os;
}

ostream& CTree::OutNode(int FromNode, int ToNode, ostream &os)	{
	int i, count;
#if TREE_DEBUG
	static vector <bool> OK;
	if(FromNode == -1) { OK.clear(); OK.assign(NoNode(),false); }
	if(OK[ToNode] == true) { cout << "\nTree error: "; OutDetail(cout); exit(-1); }
	OK[ToNode] = true;
#endif
	// If an internal node
	if(m_Node[ToNode]->m_NodeType == branch)	{
		// First visit to node gives an open parenthesis
		if(FromNode != -1) { os << "("; }	// Ensures trifurcation at root
		// Organise node visits
		count = 0; FOR(i,3)	{
			if(m_Node[ToNode]->m_viLink[i] == FromNode || m_Node[ToNode]->m_viLink[i] == -1) { continue; }
			if(count == 2) { os << ","; }	// For first node put in the extra ','
			OutNode(ToNode,m_Node[ToNode]->m_viLink[i],os);
			if(count == 0) { os << ","; }
			count++;
		}
		// The return visit gives the closing parenthesis and, if required, the branch length
		if(FromNode != -1) { os << ")"; }	// Ensures trifurcation at root
		if(FromNode != -1) { OutBranch(ToNode,FromNode,os); }
	}	else if(m_Node[ToNode]->m_NodeType == leaf)	{	// If an external node
		if(ToNode < m_iNoSeq && m_vsName.size() == m_iNoSeq) { if(m_bOutName) { os << m_vsName[ToNode]; } else { os << ToNode+1; } } else { os << ToNode + 1; }
		OutBranch(ToNode,NodeLink(ToNode,0),os);
		if(FromNode == -1) {
			if(NodeLink(ToNode,0) < m_iNoSeq && m_vsName.size() == m_iNoSeq) { if(m_bOutName) { os << "," << m_vsName[NodeLink(ToNode,0)]; } else { os << "," << NodeLink(ToNode,0); } }
			else { os << "," << NodeLink(ToNode,0) + 1; }
			if(m_iOutBra==1) { os << ":0.0"; }
		}
	} else if(m_Node[ToNode]->m_NodeType == root) {
		// First visit to node gives an open parenthesis
		if(FromNode != -1) { os << "("; }	// Ensures trifurcation at root
		// Organise node visits
		count = 0; FOR(i,2)	{
			if(m_Node[ToNode]->m_viLink[i] == FromNode || m_Node[ToNode]->m_viLink[i] == -1) { continue; }
			if(count == 2) { os << ","; }	// For first node put in the extra ','
			OutNode(ToNode,m_Node[ToNode]->m_viLink[i],os);
			if(count == 0) { os << ","; }
			count++;
		}
		// The return visit gives the closing parenthesis and, if required, the branch length
		if(FromNode != -1) { os << ")"; }	// Ensures trifurcation at root
		if(FromNode != -1) { OutBranch(ToNode,FromNode,os); }
	}
	return os;
}
ostream& CTree::OutBranch(int ToNode, int FromNode, ostream &os)	{
	int i;
	switch(m_iOutBra)	{
		case 0: break;	// No branch
		case 1:	// Simple branch length
			if(m_iNoSeq == 2 && ToNode == 1) { os << ":0.0"; break; }
			FOR(i,(int)m_Node[ToNode]->m_viLink.size()) { if(m_Node[ToNode]->m_viLink[i] == FromNode) { break; } }
			assert(i != m_Node[ToNode]->m_viLink.size());
			os << ":" << B(m_Node[ToNode]->m_viBranch[i]);
			if(m_bOutLabel && BranchLabels()) { os << "#" << m_viBranchLabels[m_Node[ToNode]->m_viBranch[i]]; }
			break;
		case 2: // Branch label
			FOR(i,(int)m_Node[ToNode]->m_viLink.size()) { if(m_Node[ToNode]->m_viLink[i] == FromNode) { break; } }
			assert(i != m_Node[ToNode]->m_viLink.size());
			os << ":b" << m_Node[ToNode]->m_viBranch[i];
			if(m_bOutLabel && BranchLabels()) { os << "#" << m_viBranchLabels[m_Node[ToNode]->m_viBranch[i]]; }
			break;
		default:
			Error("Unknown request for m_iOutBra...\n\n");
	};
	return os;
}

bool CTree::OutDetail(ostream &os, bool ForceExit)	{
	int i;
	os << "\nTree has " << m_iNoNode << " nodes and " << m_iNoBra << " branches";
	os << "\nCalculations start at Node: " << m_iStartCalc;
	os << "\nNodes: ";
	FOR(i,m_iNoNode) { os << "\n\tNode["<<i<<"]: " << *m_Node[i]; }
	os << "\nBranch: ";
	FOR(i,m_iNoBra) { os << "\n\tBranch["<<i<<"] links ("  << m_vBraLinks[i][0] << ":" << m_vBraLinks[i][1] << "): " << m_vdBra[i]; }
	if(ForceExit) { exit(-1); }
	os << "\nTree: " << *this;
	return true;
}

// Create a consistent way outputting trees
vector <int> CTree::ConstOut()	{
	int i,j;
	vector <int> List;
	FOR(i,m_iNoSeq)	{ FOR(j,i)	{ List.push_back(NodeDist(i,j)); }	}
	assert(List.size() == (((m_iNoSeq * m_iNoSeq) - m_iNoSeq) / 2));
	return List;
}

int CTree::NodeDist(int Node1, int Node2, int NodeFrom)	{
	int i;
	static int dist;
	static bool Fix;
	if(NodeFrom == -1) {
		if(Node1 == Node2) { return 0; }
		Fix = false; dist = 0;
	}
	FOR(i,(int)m_Node[Node1]->m_viLink.size())	{
		if(m_Node[Node1]->m_viLink[i] == NodeFrom || m_Node[Node1]->m_viLink[i] == -1) { continue; }
		dist++;
		if(Node2 == m_Node[Node1]->m_viLink[i]) { Fix=true; return dist; }
		NodeDist(m_Node[Node1]->m_viLink[i],Node2,Node1);
		if(Fix == true)  { return dist; }
	}
	return --dist;
}

int CTree::BranchDist(int B1, int B2, bool AllowZero)	{
	int i,j;
	int mindist = BIG_NUMBER, dist;
	if((B1 < 0 || B1 > NoBra()) || (B2 < 0 || B2 > NoBra())) { return 1; };
	FOR(i,2) {
		if(BraLink(B1,i) == -1) { continue; }
		FOR(j,2) {
			if(BraLink(B2,j) == -1) { continue; }
			dist = NodeDist(BraLink(B1,i),BraLink(B2,j));
			if(dist < mindist) { mindist = dist; }
	}	}
	if(mindist == BIG_NUMBER) { return -1; }
	return mindist;
}

void CTree::AssignNodeType(int N, ENodeType Type)	{
	assert(InRange(N,0,NoNode()));
	switch(Type)	{
	case branch:
		assert(NoLinks(N) == 3);
		m_Node[N]->m_NodeType = branch;
		break;
	case leaf:
		assert(NoLinks(N) == 1);
		m_Node[N]->m_NodeType = leaf;
		break;
	case root:
		assert(NoLinks(N) == 2);
		m_Node[N]->m_NodeType = root;
		break;
	default:;
		Error("Trying to assign a node to unknown type of ENodeType\n");
	}
}

/////////////////// Tree based pairwise distance functions //////////////////////////
// Get only the distances between the leaf nodes
vector <double> CTree::GetTreePW() {
	int i,seq;
	vector <double> dist(m_iNoSeq*m_iNoSeq,0);
	vector <double> pdist(m_iNoSeq,0);
	// Get the distances using recursive function
	FOR(seq,m_iNoSeq) {
		PWDistSub(seq,-1,&pdist);
		FOR(i,m_iNoSeq) { dist[(seq*m_iNoSeq) + i] = pdist[i]; }
	}
	return dist;
}

// Get distances between leaf and internal node distances
vector <double> CTree::GetAllTreePW()	{
	int i,j;
	bool LabelledNodes = false;
	vector <double> dists, retdist;
	dists.resize(m_iNoNode); retdist.resize(m_iNoNode*m_iNoNode);

	FOR(i,m_iNoNode) {
		// Get simple pairwise distances
		// Get node to tip distances
		if(i < m_iNoSeq) { j = i; } else {
			FOR(j,m_iNoNode) { if(m_Node[j]->m_iInternalNodeNum == i+1) { LabelledNodes = true; break; } }
			if(j == m_iNoNode) { j = i; if(LabelledNodes == true) { cout << "\nWarning: Tree::GetAllTreePW() some nodes labelled whilst some are not..."; } }
		}
		PWDistSub(j,-1,&dists,true);
		FOR(j,m_iNoNode) { retdist[(i*m_iNoNode)+j] = dists[j]; }
	}
//	cout << "\nThe node calculations are:"; FOR(i,m_iNoNode) { cout << "\nNode["<<i<<"]:"; FOR(j,m_iNoNode) { cout << "\t" << retdist[(i*m_iNoNode)+j]; } }
	return  retdist;
}

void CTree::PWDistSub(int NodeTo, int NodeFrom,vector <double> *d, bool DoInternalNodes) {
	int i;
	static double dist;
	// Initialise if first node
	if(NodeFrom == -1) {
		dist = 0.0;
		assert((DoInternalNodes == false && d->size() == m_iNoSeq) || (DoInternalNodes == true && d->size() == m_iNoNode));
		FOR(i,(int)d->size()) { d->at(i) = -1; }
	}
	// If a terminal node then store the distance
	if(NodeTo < m_iNoSeq) { d->at(NodeTo) = dist; } else if(DoInternalNodes == true) { d->at(NodeTo) = dist; }
	// recurse down tree
	FOR(i,(int)m_Node[NodeTo]->m_viLink.size()) {
		if(m_Node[NodeTo]->m_viLink[i] == NodeFrom || m_Node[NodeTo]->m_viLink[i] == -1) { continue; }
		dist += B(m_Node[NodeTo]->m_viBranch[i]);
		PWDistSub(m_Node[NodeTo]->m_viLink[i],NodeTo,d,DoInternalNodes);
		dist -= B(m_Node[NodeTo]->m_viBranch[i]);
	}
}

// Recursive function that collects the nodes of exactly NodeDepth, will also include leaf nodes if less than NodeDepth if GetLess == true
vector <int> CTree::GetNodesOfDepth(int InitNode, int NodeDepth, bool GetLess, vector <int> *NodesFrom, vector <int> *NodeCov, vector <double> *ExtBra,int NodeFr, bool First)	{
	static vector <int> vNodes;
	static int CurDepth;
	int i;
	// Check entry conditions
	assert(InRange(InitNode,0,m_iNoNode) && InRange(NodeFr,-1,m_iNoNode) && (NodeDepth > 0));
	// Initialise if required
	if(First == true) { vNodes.clear(); CurDepth = 0; if(NodesFrom != NULL) { NodesFrom->clear(); } if( NodeCov != NULL) { NodeCov->clear(); } }
	// Get the nodes if the node depth is reached
	if(CurDepth == NodeDepth)	{
		vNodes.push_back(InitNode);
		if(NodesFrom != NULL) { NodesFrom->push_back(NodeFr); }
		// Find Branch
		if(ExtBra != NULL)	{
			FOR(i,m_iNoBra) {
				if( (m_vBraLinks[i][0] == InitNode && m_vBraLinks[i][1] == NodeFr) || (m_vBraLinks[i][0] == NodeFr && m_vBraLinks[i][1] == InitNode) )	{
					ExtBra->push_back(m_vdBra[i]); break;
			}	}
			if(i == m_iNoBra) { Error("\nCannot find branch in GetNodesOfDepth...\n"); }
		}
		return vNodes;
	}
	// Get the covered nodes
	if(NodeCov != NULL) {
		if(!IsIn(InitNode,*NodeCov) && InitNode >= m_iNoSeq) { NodeCov->push_back(InitNode); }
	}
	// Check for leaf nodes
	if(InitNode < m_iNoSeq)		{
		if(GetLess == true) {
			vNodes.push_back(InitNode);
			if(NodesFrom != NULL) { NodesFrom->push_back(NodeFr); }
			// Find Branch
			if(ExtBra != NULL)	{
				FOR(i,m_iNoBra) {
					if( (m_vBraLinks[i][0] == InitNode && m_vBraLinks[i][1] == NodeFr) || (m_vBraLinks[i][0] == NodeFr && m_vBraLinks[i][1] == InitNode) )	{
						ExtBra->push_back(m_vdBra[i]); break;
				}	}
				if(i == m_iNoBra) { Error("\nCannot find branch in GetNodesOfDepth...\n"); }
		}	}
		return vNodes;
	}
	// Otherwise descend the tree by looping through the nodes
	CurDepth++;
	FOR(i,(int)m_Node[InitNode]->m_viLink.size()) {
		if(m_Node[InitNode]->m_viLink[i] == NodeFr || m_Node[InitNode]->m_viLink[i] == -1) { continue; }
		GetNodesOfDepth(m_Node[InitNode]->m_viLink[i],NodeDepth,GetLess,NodesFrom,NodeCov,ExtBra,InitNode,false);
	}
	CurDepth--;
	return vNodes;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Centre-Point algorithms for SNAP

// Branch centre point
vector <int> CTree::BranchCP(int CP, int Depth, vector <int> *NodeRem, vector <int> *NodeCovered, vector <double> *ExtBra)	{
	int i,j;
	vector <int> Nodes, temp, RetNode, NC;
	NodeRem->clear(); if(NodeCovered != NULL) { NodeCovered->clear(); }
	// Clean in preparation of new subtree
	FOR(i,2) {
		j = 0; if(i==0) { j = 1; }
		Nodes = GetNodesOfDepth(m_vBraLinks[CP][i],Depth,true,&temp,&NC,ExtBra,m_vBraLinks[CP][j]);
		assert(Nodes.size() == temp.size());
		// Store the nodes for the leafmap and the nodes from for calculating partial likelihoods
		FOR(j,(int)Nodes.size()) { RetNode.push_back(Nodes[j]); NodeRem->push_back(temp[j]); }
		if(NodeCovered != NULL) { FOR(j,(int)NC.size()) { if(!IsIn(NC[i],*NodeCovered)) { NodeCovered->push_back(NC[j]); } } }
	}
	return RetNode;
}
// Node centre point
vector <int> CTree::NodeCP(int Node, int Depth, vector <int> *NodeRem, vector <int> *NodeCovered, vector <double> *ExtBra)	{
	vector <int> Nodes, NodeFr, NC;
	// Clean the node lists
	NodeRem->clear(); if(NodeCovered != NULL) { NodeCovered->clear(); }
	// Clean in preparation of new subtree
	Nodes = GetNodesOfDepth(Node,Depth,true,NodeRem,NodeCovered,ExtBra);
	assert(Nodes.size() == NodeRem->size());
	return Nodes;
}

// Tree replacement function
void CTree::ReplaceTreeCP(CTree *NT,vector <int> LeafMap,vector <int> NCover, bool VerifyBranchLinks)	{
	int i,j,nf, nt;
	vector <int> New2Old, Bran, NewL, NewB, NodeIntVals;
	assert(LeafMap.size() + NCover.size() == NT->m_iNoNode);
	New2Old = LeafMap;
#if TREE_DEBUG
	int count
	// Check entry conditions
	FOR(i,LeafMap.size())	{ // Check all leaves map to one and only one NCover
		count = 0;
		FOR(j,NoLinks(LeafMap[i])) { if(IsIn(NodeLink(LeafMap[i],j),NCover)) { count ++; } }
		assert(count == 1);
		if(count != 1) { cout << "\nCTree::ReplaceTreeCP -- Error in LeafList " << LeafMap[i] << "; count = " << count << "..."; OutDetail(); exit(-1); }
	}
	FOR(i,NCover.size())	{ // Check all nodes covered link only to other nodes covered and leafmaps
		count = 0;
		FOR(j,NoLinks(NCover[i])) { if(IsIn(NodeLink(NCover[i],j),NCover) || IsIn(NodeLink(NCover[i],j),LeafMap)) { count ++; } }
		assert(count == NoLinks(NCover[i]));
		if(count != NoLinks(NCover[i])) { cout << "\nCTree::ReplaceTreeCP -- Error in NCover " << NCover[i] << "..."; OutDetail(); exit(-1); }
	}
#endif
	FOR(i,(int)NCover.size()) {
		// Get branches
		FOR(j,(int)m_Node[NCover[i]]->m_viBranch.size()) {
			if(!IsIn(m_Node[NCover[i]]->m_viBranch[j],Bran)) { Bran.push_back(m_Node[NCover[i]]->m_viBranch[j]); }
		}
		// Get Node nums and delete them
		New2Old.push_back(NCover[i]); NodeIntVals.push_back(m_Node[NCover[i]]->m_iInternalNodeNum);
		delete m_Node[NCover[i]]; m_Node[NCover[i]] = NULL;
	}
	sort(Bran.begin(),Bran.end());
	assert(Bran.size() == NT->m_iNoBra);
	// Apply the nodes
	FOR(nf,NT->m_iNoNode) {
		nt = New2Old[nf];		// Set NodeTo;
		// Do leaf nodes
		if(nf<NT->m_iNoSeq) {
			// find the link that goes towards another covered node
			FOR(i,(int)m_Node[nt]->m_viLink.size()) { if(IsIn(m_Node[nt]->m_viLink[i],New2Old)) { break; } }
			assert(i != m_Node[nt]->m_viLink.size() && NT->m_Node[nf]->m_viLink.size() == 1);
			// Set the link
			m_Node[nt]->m_viLink[i] = New2Old[NT->m_Node[nf]->m_viLink[0]];
			// Set the branch
			m_Node[nt]->m_viBranch[i] = Bran[NT->m_Node[nf]->m_viBranch[0]];
		} else {		// Do internal nodes
			NewL.clear(); NewB.clear();
			// Sort out the links
			FOR(i,(int)NT->m_Node[nf]->m_viLink.size()) {
				NewL.push_back(New2Old[NT->m_Node[nf]->m_viLink[i]]);
				NewB.push_back(Bran[NT->m_Node[nf]->m_viBranch[i]]);
			}
			// Create the new node
			assert(InRange( (int) (nf - NT->m_iNoSeq) ,(int) 0,(int) NodeIntVals.size()));
			if(m_Node[nt] != NULL) { delete m_Node[nt]; }
			m_Node[nt] = new CNode(NewL,NewB,NodeIntVals[nf - NT->m_iNoSeq]);
			// Apply the branch lengths
			FOR(i,(int)m_Node[nt]->m_viBranch.size()) {
				SetB(m_Node[nt]->m_viBranch[i],NT->B(NT->m_Node[nf]->m_viBranch[i]),true,true);
	}	}	}
	BuildBraLinks(VerifyBranchLinks);
	OrderNode();
	m_bFastCalcOK = false;
	ValidateTree();
}


///////////////////////////////////////////////////////////////////////////////
// Gets the base pairwise distances for a tree without the subtree (defined by LeafMap & NCover)
// branches counted.
vector <double> CTree::GetPartialTreeDist(vector <int> LeafMap, vector <vector <int> > NBelow)	{
	int Leaf,Leaf2,i,j;
	vector <double> TreeDist(m_iNoSeq*m_iNoSeq,0);
	vector <double> dist(m_iNoSeq,0);
	// Check entry conditons
	assert(NBelow.size() == LeafMap.size());
	// Do the calculations
	FOR(Leaf,(int)LeafMap.size())	{	// Loop through the leaves and get the pairwise distances
		if(NBelow[Leaf].size() == 1) { continue; }
		PWDistSub(LeafMap[Leaf],-1,&dist);	// Distances between LeafMap nodes and the true leaves
		assert(dist.size() == m_iNoSeq);
		// Loop through the Nodes Below and fill in the partial distances
		FOR(i,(int)NBelow[Leaf].size()) {	// This node
			FOR(Leaf2,(int)LeafMap.size()) {
				if( (NBelow[Leaf].size() == 1 && NBelow[Leaf2].size() == 1) || Leaf==Leaf2) { continue; }
				FOR(j,(int)NBelow[Leaf2].size())	{
					TreeDist[(m_iNoSeq*NBelow[Leaf][i]) + NBelow[Leaf2][j]] += dist[NBelow[Leaf][i]];
					TreeDist[(m_iNoSeq*NBelow[Leaf2][j]) + NBelow[Leaf][i]] += dist[NBelow[Leaf][i]];
		}	}	}
		// Now do the other distances that don't span the subtree
		FOR(i,(int)NBelow[Leaf].size())	{
			PWDistSub(NBelow[Leaf][i],-1,&dist);
			FOR(j,i)	{
				TreeDist[(m_iNoSeq*NBelow[Leaf][i])+NBelow[Leaf][j]] += dist[NBelow[Leaf][j]];
				TreeDist[(m_iNoSeq*NBelow[Leaf][j])+NBelow[Leaf][i]] += dist[NBelow[Leaf][j]];
	}	}	}
	return TreeDist;
}

//////////////////////////////////////////////////////////////
// Get full pairwise distances from a subtree and a set of
// partial RMSDs taken from GetPartialTreeDist(...) above
vector <double> CTree::GetSubTreePW(vector <int> LeafMap, vector <vector <int> > NBelow, vector <double> Dist)	{
	int Leaf,Leaf2,i,j,OriNoSeq = 0;
	vector <double> SubDist = GetTreePW();	// The pairwise distances for the subtree
	// Count the number of sequences
	FOR(i,(int)NBelow.size()) { OriNoSeq += (int)NBelow[i].size(); }
	// Check entry conditions
	assert(LeafMap.size() == NBelow.size()); assert(Dist.size() == OriNoSeq * OriNoSeq);
	// Do the distances
	FOR(Leaf,(int)LeafMap.size())	{
		FOR(i,(int)NBelow[Leaf].size())	{
			FOR(Leaf2,Leaf)	{
				FOR(j,(int)NBelow[Leaf2].size()) {
					Dist[(NBelow[Leaf][i] * OriNoSeq)+NBelow[Leaf2][j]] += SubDist[(Leaf*m_iNoSeq)+Leaf2];
					Dist[(NBelow[Leaf2][j] * OriNoSeq)+NBelow[Leaf][i]] += SubDist[(Leaf*m_iNoSeq)+Leaf2];
	}	}	}	}
	return Dist;
}

int CTree::BestStartCalc()	{
	int i, BestBra = -1;
	double temp = -1.0, BScore = -1.0;
	// Find the most efficient place to perform the calculation
	FOR(i,NoBra()) { 	temp = QuadB(i); if(temp > BScore) { BScore = temp; BestBra = i; } }
	assert(InRange(BestBra,0,NoBra()));
	m_iStartCalc = BestBra;
	m_bFastCalcOK = false;
	return m_iStartCalc;
}

///////////////////////////////////////////////////////////////////////
// Routines to build subtrees

void CTree::BuildOriSubTree(CTree *T, vector <bool> NodesBool)	{
	int i;
	vector <int> NC;						// Values for passing to another BuildOriSubTree routine
	assert(NodesBool.size() == NoNode() + NoSeq());
	FOR(i,(int)NodesBool.size())	{ if(NodesBool[i] == true) { NC.push_back(i); } }
	BuildOriSubTree(T,NC);
}

void CTree::BuildOriSubTree(CTree *T, vector <int> NC)	{
	int i,j, NodCount, NPoint;
	vector <int> NCover,LeafMap,NFrom;			// The values that NC shall be changed to
	// Break up NC into the correct numbers
	FOR(i,(int)NC.size())	{
		if(NC[i] < NoSeq()) { LeafMap.push_back(NC[i]); NFrom.push_back(NodeLink(NC[i],0)); continue; }	// Do leaf nodes
		// Find out whether the node is a leaf node or an internal node
		NodCount = 0; NPoint = -1;
		FOR(j,NoLinks(NC[i])) {
			if(IsIn(NodeLink(NC[i],j),NC)) {
				NodCount++; NPoint = NodeLink(NC[i],j);
				if(NodCount > 1) { break; }
		}	}
		// NodCount == 1 means that it is a leaf node
		if(NodCount == 1) { LeafMap.push_back(NC[i]); NFrom.push_back(NPoint); }
		// NodCount > 1 means its a covered node
		else if(NodCount > 1) { NCover.push_back(NC[i]); }
		// NodCount == 0 is an error
		else { Error("Error when constructing subtree in CTree::BuildOriSubTree(...)\n"); }
	}
	assert(LeafMap.size() == NFrom.size());
	BuildOriSubTree(T,LeafMap,NCover,NFrom);
}

void CTree::BuildOriSubTree(CTree *RetTree, vector <int> LeafMap, vector <int> NCover, vector <int> NFrom)	{
	int i,j,k,Node,NoSeq = (int)LeafMap.size();
	vector <int> Link, Branch; Link.assign(3,-1); Branch.assign(3,-1);
	vector <double> BVal; BVal.assign(3,-1.0);
	bool flag;
	assert(!LeafMap.empty() && LeafMap.size() == NCover.size() + 2);
	// Do only Nodes
	RetTree->GetMemory(NoSeq);
	FOR(i,NoSeq)	{	// LeafNodes (easy)
		// Find the parent node
		FOR(Node,(int)NCover.size()) { if(NCover[Node] == NFrom[i]) { break; } }
		assert(Node != NCover.size());
		// Create the child node
		RetTree->m_Node[i]->SetNode(Node + NoSeq,i,-1);
	}
	FOR(i,(int)NCover.size())	{	// Branch nodes
		FOR(j,3)	{	// Get each of the three links
			// Find the real tree node
			Node = m_Node[NCover[i]]->m_viLink[j];
			// Map back to the subtree nodes
			flag = false; FOR(k,NoSeq) { if(LeafMap[k] == Node) { flag = true; break; } }
			if(k == NoSeq) {
				FOR(k,(int)NCover.size()) { if(NCover[k] == Node) { k += NoSeq; flag = true; break; } }
			}
			assert(flag == true);
			Link[j] = k;
		}
		sort(Link.begin(),Link.end());
		RetTree->m_Node[i+NoSeq]->SetNode(Link[0],Link[1],Link[2],-1,-1,-1,-1);
	}
	// Now build the branches and their links
	RetTree->BuildBranches();
	RetTree->BuildBraLinks();
	// Now copy branches lengths from original tree to RetTree
	Link.clear(); Link = LeafMap; FOR(i,(int)NCover.size()) { Link.push_back(NCover[i]); }
	FOR(i,(int)RetTree->m_vdBra.size())	{
		RetTree->SetB(i,B(FindBra(Link[RetTree->m_vBraLinks[i][0]],Link[RetTree->m_vBraLinks[i][1]])),true);
	}
	RetTree->m_bReady = true;
}

//////////////////////////////////////////////////////////////////////////
// Function used for investigating knots in the tree
vector <vector <int> > CTree::GetKnotClusters(vector <bool> INC, int ChangeRad)	{
	int i,j,k;
	vector <vector <int> > Clusters;
	vector <int> v_C,v_Blank;
	assert(INC.size() == NoNode());

	// Get original clusters
	for(i=NoSeq();i<NoNode();i++)	{
		if(INC[i] == false) { continue; }
		v_Blank = GetNodesOfDepth(i,ChangeRad,true,NULL,&v_C,NULL);
		FOR(j,(int)v_Blank.size()) { v_C.push_back(v_Blank[j]); }
//			Nodes = GetNodesOfDepth(m_ariBraLinks[CP][i],Depth,true,&temp,&NC,ExtBra,m_ariBraLinks[CP][j]);
		// See whether the nodes need a new cluster
		FOR(j,(int)Clusters.size()) {
			FOR(k,(int)v_C.size()) {
				if(IsIn(v_C[k],Clusters[j])) { // Is in cluster so store
					FOR(k,(int)v_C.size()) { if(!IsIn(v_C[k],Clusters[j])) { Clusters[j].push_back(v_C[k]); } }
					k = -1; break;
			}	}
			if(k == -1) { break; }
		}
		if(k != -1) { Clusters.push_back(v_C); }
		v_C.clear();
	}
//	cout << "\nFinished with " << Clusters.size() << " clusters";
//	FOR(i,Clusters.size()) { Sort(&Clusters[i]); cout << "\nCluster["<<i<<"] size = " << Clusters[i].size()<< ": " << Clusters[i]; }
//	exit(-1);
	return Clusters;
}

//////////////// Function applying names to a tree ///////////////////
void CTree::SetNames(vector <string > NewNames,bool Overwrite) {
	if(!Names().empty() && Overwrite == false) { Error("\nTrying to overwrite names when not permitted in CTree::SetNames(...)\n"); }
	if(NewNames.size() != m_iNoSeq) { Error("\nTrying to write " + int_to_string(NewNames.size()) + " NewNames to a tree with " + int_to_string(m_iNoSeq) + " sequences in CTree::SetNames\n"); }
	m_vsName = NewNames;
}



bool IsSameTree(CTree *T1, CTree *T2)	{
	assert(!T1->IsCutTree() && !T2->IsCutTree());
	int i,j;
	vector < vector <int> > L1;
	vector <int> L2,R2;
	// Do the obvious checks first!
	if(T1->NoSeq() != T2->NoSeq()) { return false; }
	// Get branch sets for T1
	FOR(i,T1->NoBra())	{
		L2.clear(); R2.clear();
		T1->BranchSets(i,L2,R2);
		L1.push_back(L2);
	}
	// Now compare these branch sets to those of T2
	FOR(i,T2->NoBra()) {
		L2.clear(); R2.clear(); T2->BranchSets(i,L2,R2); // Get branch sets for Tree 2
		// See if they match T1
		FOR(j,T1->NoBra()) { if(Compare(&L1[j],&L2)) { break; } }
		if(j == T1->NoBra()) { return false; }
	}
	// If these conditions are all met, then they're the same tree
	return true;
}

/////////////////////////////////////////////////////////////////////////////////
// Functions for finding maximum tree length subtree containing exactly SubSeq
// from the tree FullTree
CTree FindGreedySubTree(CTree *FullTree, int SubSeq) {
	int i,j,max,bra, best_seq, best_bra;
	double dist, best_val;
	vector <double> PWDists;
	vector <int> SeqsToAdd; FOR(i,FullTree->NoSeq()) { SeqsToAdd.push_back(i); }
	CTree CurTree;
	string Start;
	// Get starting pair of sequences
	// 0. Need to ensure there are no real multifurcations because these can end up being resolved incorrectly.
	FOR(i,FullTree->NoBra()) {
		if(FullTree->B(i) < FLT_EPSILON) { FullTree->SetB(i,Random()*1.0E-4); }
	}
	// 1. Get PW distances
	PWDists = FullTree->GetTreePW();
	// 2. Find the minimum pairwise distance and initialise the tree
	dist = 0; max = 0; FOR(i,(int)PWDists.size()) { if(PWDists[i] > dist) { dist = PWDists[i]; max = i; } }
	Start = "(" + int_to_string((max / FullTree->NoSeq())+1) + ":" + double_to_string(PWDists[max]) + "," + int_to_string((max % FullTree->NoSeq())+1) + ":0.0);";
	SeqsToAdd[max / FullTree->NoSeq()] = -1; SeqsToAdd[max % FullTree->NoSeq()] = -1;
	CurTree.CreateTree(Start,FullTree->Names(),true,true,true);
	// Progressively add the sequences to the tree
//	FullTree->OutBra(); cout << "\nInitial starting tree:\n" << *FullTree;
	PWDists = FullTree->GetAllTreePW();
	FOR(i,SubSeq-2)	{
//		cout << "\nDoing sequence add #" << i << "\n\tTrying";
		// Find which sequence gives greatest gain
		best_seq = -1; best_val = -BIG_NUMBER;
		FOR(j,(int)SeqsToAdd.size()) {
//			cout << " " << SeqsToAdd[j] << flush;
			if(SeqsToAdd[j]<0) { continue; }
			// Initialise the tree part
			dist = TravAddGreedy(&CurTree,(int)CurTree.StartCalc(),-1,SeqsToAdd[j],&PWDists,&bra);
//			cout << "^" << flush;
			if(dist > best_val) { best_seq = j; best_val = dist; best_bra = bra; }
		}
//		cout << "\n\tBest " << SeqsToAdd[best_seq] << flush;
		GreedySeq2Tree(best_bra,SeqsToAdd[best_seq],&CurTree,&PWDists); SeqsToAdd[best_seq] = -1;
//		CurTree.OutBra(); cout << "\n\t" << CurTree;
		// Check compatibility between new and full tree
		assert(FullTree->IsCompatible(CurTree));
	}
	return CurTree;
}

// Core in order tree traversal for adding a sequence
// Need to find the closest location to add a sequence
double TravAddGreedy(CTree *CurT, int To, int Fr, int Seq2Add, vector <double> *PWDist, int *BestBra) {
	int i;
	static double best;	// Current best value
	double dist;
	if(CurT->NodeType(To) == leaf) {	// Leaf nodes
		if(Fr == -1) { // Starting node
			*BestBra = CurT->FindBra(To,CurT->NodeLink(To,0));
			best = GetDist(Seq2Add,To,CurT->NodeLink(To,0),PWDist); // Always store the first as best
		} else {
			dist = GetDist(Seq2Add,To,Fr,PWDist);
			if(dist < best) { *BestBra = CurT->FindBra(To,Fr); best = dist; }
		}
	} else { 						// internal nodes
		FOR(i,CurT->NoLinks(To))	{ if(CurT->NodeLink(To,i) == Fr || CurT->NodeLink(To,i) == -1) { break; } }
		assert(i != CurT->NoLinks(To));
		// If the node from isn't a leaf node do the internal calculation (i.e. avoids first node)
		if(CurT->NodeType(CurT->NodeLink(To,i)) != leaf)	{
			dist = GetDist(Seq2Add,To,Fr,PWDist);
			if(dist < best) { *BestBra = CurT->FindBra(To,Fr); best = dist; }
	}	}
	// Do the looping
	FOR(i,CurT->NoLinks(To))	{
		if(CurT->NodeLink(To,i) == Fr || CurT->NodeLink(To,i) == -1) { continue; }
		TravAddGreedy(CurT,CurT->NodeLink(To,i),To,Seq2Add, PWDist,BestBra);
	}
	return best;
}
// Get the distance obtained from adding a sequence (Add) between nodes (a,b)
double GetDist(int Add, int a, int b, vector <double> *PWDist) {
	int index = sqrt(PWDist->size());
	double d_ab = PWDist->at((a*index)+b), d_aAdd = PWDist->at((a*index)+Add),d_bAdd = PWDist->at((Add*index)+b), l = d_ab+d_aAdd+d_bAdd;
	return 0.5 * (l - (2*d_ab));

}

void GreedySeq2Tree(int Br, int Seq2Add, CTree *CurTree, vector <double> *PWDist) {
	int i,j,k;
	int tmp, NLeft = CurTree->BraLink(Br,0), NRight = CurTree->BraLink(Br,1);
	vector <int> vtmp;
	if(NRight > NLeft) { tmp = NRight; NRight = NLeft; NLeft = tmp; }
	double Dleft = GetDist(NLeft,NRight,Seq2Add,PWDist), Dright = GetDist(NRight,NLeft,Seq2Add,PWDist), Dnew = GetDist(Seq2Add,NLeft,NRight,PWDist);
	double prop;
	// Error check and get the branch proportions
//	cout << "\nDleft (" << Dleft << ") + Dright ("<< Dright << ") == " << Dleft + Dright << " cf. " << CurTree->B(Br) << " diff: " << fabs(CurTree->B(Br) - (Dleft + Dright));
	assert(fabs((Dleft + Dright) - CurTree->B(Br)) < 10-4);	// Check branch lengths agree
	prop = Dright / (Dleft + Dright);
	// Create the new node and branch for adding to the tree
	// 1. Find a spare node and branch
	assert(CurTree->NoLinks(Seq2Add) == 0);	// Check the node is empty
	// Find the node that matches the expected distances
	// Note this this not guarantee the correct mapping between the CurTree and original tree due to multifurcations, but it works!
	for(i=CurTree->NoSeq();i<CurTree->NoNode();i++) {	// Find spare internal node
		if(CurTree->NoLinks(i) == 0)  {
			// Check 3 way distance condition matches
			if(fabs(PWDist->at((CurTree->NoNode() * i) + NLeft) - Dleft) < FLT_EPSILON &&
				fabs(PWDist->at((CurTree->NoNode() * i) + NRight) - Dright) < FLT_EPSILON &&
				fabs(PWDist->at((CurTree->NoNode() * i) + Seq2Add) - Dnew) < FLT_EPSILON) { break; }
	} 	}
	assert(i != (int)CurTree->NoNode());
	// Find two spare branches
	k = -1; FOR(j,CurTree->NoBra()) { if(!CurTree->GoodBra(j)) { if(k==-1) { k = j; } else { /*cout << "\nBranch[" <<j << "] is spare";*/ break; } } }
	// 2. Now connect them up (i = node; j = branch; k = branch spare)
	vtmp.clear(); vtmp.push_back(i); CurTree->ReplaceNodeLink(Seq2Add,vtmp); 											// } replace nodes
	vtmp.clear(); vtmp.push_back(Seq2Add); vtmp.push_back(-1); vtmp.push_back(-1); CurTree->ReplaceNodeLink(i,vtmp);	// }
	CurTree->AssignNodeType(Seq2Add,leaf); CurTree->AssignNodeType(i,branch);											// }
	vtmp.clear(); vtmp.push_back(j); CurTree->ReplaceNodeBra(Seq2Add,vtmp);		// }
	vtmp.push_back(k); vtmp.push_back(-1); CurTree->ReplaceNodeBra(i,vtmp);		// } Do branches
	CurTree->ReplaceBraLink(j,0,min(Seq2Add,i)); CurTree->ReplaceBraLink(j,1,max(Seq2Add,i));	// } Do links
	CurTree->ReplaceBraLink(k,0,i); CurTree->SetB(j,Dnew);																// }

	CurTree->AddSeq(Seq2Add,Br,prop);
}

vector <string> ReadTreeNames(string Tree) {
	vector <string> retNames;
	for(int i = 0; i < Tree.size(); i++) {
		if(Tree[i] == '(' || Tree[i] == ',') {
			if(Tree[i+1] == '(') { continue; }
			i++;
			for(int j = i; j < Tree.size(); j++) {
				if(Tree[j] == ',' || Tree[j] == ')' || Tree[j] == ':') { retNames.push_back(Tree.substr(i,j-i)); break; }
			}
		}
	}
	return retNames;
}

////////////////// Tools /////////////
bool FlipBool(bool V)   { if(V == true) { return false; } return true; }
bool FlipBin(int i)             { assert(i==0||i==1); if(i==0) { return 1; } return 0; }
string int_to_string(int num) { stringstream ss; ss << num; return ss.str(); }
string double_to_string(double num) { stringstream ss; ss << num; return ss.str(); }
ostream& Error(string str , ostream  &os) { os << str; assert(0); exit(-1); };
