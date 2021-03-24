#include "head.h"

string loc1 = "D:\\cit-HepTh\\Cit-HepTh";
string loc2 = ".txt";
string loc = loc1 + loc2;
string loc_fixed = loc1 + "-fixed" + loc2;
string loc_dic = loc1 + "-dic" + loc2;
static int num_xiaozhi = 0;
static int num_max_outEdges = 0;

typedef property<vertex_name_t, int> VertexProperties;
typedef property<edge_finished_t, int> EdgeProperties;
typedef adjacency_list<vecS, vecS, directedS, VertexProperties, EdgeProperties> Graph;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::vertex_descriptor Vertex;



class Tree;
class TreeNode  //表示树的结点，左孩子右兄弟
{
	friend Tree;
public:
	int data;//数据
	TreeNode* left;
	TreeNode* right;  //子女和兄弟指针
	TreeNode(int value = 0, TreeNode* l = NULL, TreeNode* r = NULL) :data(value), left(l), right(r) {}//构造函数
};
class Tree
{
public:
	Tree() { root = current = NULL; }; //建立树
	void setRoot(int value);
	TreeNode* getRoot();
	int toLeft();
	int toRight();
	int Parent();
	int Find(int target);
	TreeNode* root, * current;    //表示根和当前指针
//void PreOrder(ostream& out, TreeNode* p)
	int Find(TreeNode* p, int target);
	//void RemovesubTree(TreeNode* p)
	int FindParent(TreeNode* t, TreeNode* p);
};

void Tree::setRoot(int value)
{
	root = current = new TreeNode(value);  //建立根结点
}

TreeNode* Tree::getRoot()  //用来查找根,使之成为当前结点
{
	if (root == NULL) { current = NULL; return NULL; }
	else { current = root; return root; }
}

int Tree::toLeft()
{ //寻找当前结点的第一个子女.若无子女,则函数返回0,当前指针为NULL;
  //否则返回1,当前指针移到当前结点的第一个子女;
	if (current != NULL && current->left != NULL)
	{
		current = current->left; return 1;
	}
	current = NULL; return 0;
}


int Tree::toRight()
{
	//寻找当前结点的下一个兄弟.若无下一个兄弟,则函数返回0,当前指针为NULL;否则返回1,当前指针移动到当前结点的下一个兄弟.
	if (current != NULL && current->right != NULL)
	{
		current = current->right; return 1;
	}
	current = NULL; return 0;
}

int Tree::Parent()
{//寻找当前结点的双亲结点.若无双亲结点,则函数返回0,当前指针为NULL;否则返回1,当前指针移到当前结点的双亲结点
	TreeNode* p = current;//p保存当前指针
	TreeNode* t;
	if (current == NULL || current == root) { current = NULL; return 0; }//没有双亲
	t = root;   //t从根结点开始
	int  k = FindParent(t, p);  //从根结点t后序遍历找结点p的双亲结点
	return k;
}

//t传进来就是指向root
int Tree::FindParent(TreeNode* t, TreeNode* q)
{//私有函数:从根指针t所指结点后序遍历,找结点p的双亲结点,记在 current中 
	TreeNode* p = t->left;//找第一棵子树
	while (q != NULL && q != p)       //q==NULL,无子树;q==p,找到双亲
	{
		int i = FindParent(q, p);
		if (i != 0) return i;
		q = q->right; //找下一棵子树
	}
	if (q != NULL && q == p) { current = t; return 1; } //成功返回1
	else return 0;    //否则返回0
}

int Tree::Find(int target)
{//在树中搜索符合给定值target的结点.搜索成功,函数返回1,否则返回0
	if (root == NULL) return 0;      //为空树
	return Find(root, target);    //调用私有函数Find搜索
}

int Tree::Find(TreeNode* p, int target)
{ //私有函数:在由根指针p所指的树或子数中搜索其数据成员等于target的结点.搜索成功,函数返回1,同时该结点成为当前结点;否则函数返回0
	int result = 0;
	if (p->data == target) { result = 1; current = p; } //搜索成功,p成为当前结点
	else           //否则,继续搜索
	{
		TreeNode* q = p->left;
		while (q != NULL && !(result = Find(q, target))) q = q->right;
	}
	return result;
}




//保存字典
void preProcess_Cit_Hepth() {
	ifstream infile;
	//加入顶点和边
	infile.open(loc, ios::in);
	ofstream outfile_fixed;
	ofstream outfile_dic;
	outfile_fixed.open(loc_fixed,ios::out|ios::trunc);
	outfile_dic.open(loc_dic, ios::out | ios::trunc);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	string u, v = "";
	//从第五行开始取数
	int begin = 5;
	string str;
	map<string, int> dic;
	int index_dic = 0;
	map<string, int>::iterator finder;

	while (getline(infile,str)) {            // 若未到文件结束一直循环
		if (begin > 1) {
			begin--;
			continue;
		}
		istringstream strstream(str);
		strstream >> u >> v;

		finder = dic.find(u);
		if (finder == dic.end()) {
			dic.insert(pair<string, int>(u, index_dic));
			outfile_dic << u << '\t'<<index_dic<<endl;
			index_dic++;
		}

		finder = dic.find(v);
		if (finder == dic.end()) {
			dic.insert(pair<string, int>(v, index_dic));
			outfile_dic << v << '\t' << index_dic << endl;
			index_dic++;
		}
		outfile_fixed << dic[u]<< '\t' << dic[v]<< endl;

	}
	infile.close();   //关闭文件
	outfile_fixed.close();
	outfile_dic.close();

}

void in_degree_zero_DFS(int u, Graph* g, Tree *tree,TreeNode* currentNode) {
	Graph::out_edge_iterator outedgeIt, outedgeEnd;
	tie(outedgeIt, outedgeEnd) = out_edges(u, *g);
	if (!tree->getRoot()) {
		tree->setRoot(u);
		currentNode = tree->getRoot();
		cout << currentNode->data << endl;
	}
	for (; outedgeIt != outedgeEnd; ++outedgeIt) {
		Edge e = *outedgeIt;
		property_map<Graph, edge_finished_t>::type edge_finished1 = get(edge_finished, *g);
		Vertex v_tar = target(e, *g);
		if (((Tree)*tree).Find(v_tar)| edge_finished1[e])
			continue;
		TreeNode* treeNode = new TreeNode(v_tar);
		//给currentNode->left加兄弟，currentNode->left是不是NULL都无所谓
		TreeNode* temp = currentNode->left;
		currentNode->left = treeNode;
		treeNode->right = temp;
		//g中标记这条边
		edge_finished1[e] = 1;
		in_degree_zero_DFS(v_tar, g, tree, currentNode);
	}

}

void graph_proc() {
	
	//G=(V,E,L_V)
	//T={V,E,S_T(T中顶点的数量),L_V}
	//使用bidirectionalS 替换directedS，这样可以使用in_edge、in_degree等，但消耗要增一倍
	Graph g;
	int x, y = 0;
	string str;
	ifstream infile;
	//加入顶点和边
	infile.open(loc_fixed, ios::in);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {            // 若未到文件结束一直循环
		istringstream strstream(str);
		strstream >> x >> y;
		add_edge(x, y, g);
	}
	infile.close();   //关闭文件
	int num_v = num_vertices(g);
	int* dfn= new int[num_v];
	int* low= new int[num_v];
	bool* vis= new bool[num_v];
	int* stack= new int[num_v];
	Tree** arr_xiaozhi =new Tree *[num_v];

	int* in_degree = new int[num_v];
	for (int i = 0; i < num_v; i++)
		in_degree[i] = 0;
	//加入顶点和边
	infile.open(loc_fixed, ios::in);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {            // 若未到文件结束一直循环
		istringstream strstream(str);
		strstream >> x >> y;
		in_degree[y] = in_degree[y] + 1;
	}

	infile.close();   //关闭文件
	//Graph::out_edge_iterator outedgeIt, outedgeEnd;
	//tie(outedgeIt, outedgeEnd) = out_edges(0, g);
	//for (; outedgeIt != outedgeEnd; ++outedgeIt)
	//{
	//	Edge e = *outedgeIt;
	//	Vertex v = target(e, g);
	//}

	//遍历所有顶点，将入度为0顶点的全部标记
	//用bgl太慢了，直接建数组存
	for (int i = 0; i < num_v; i++) {
		num_max_outEdges > in_degree[i] ? NULL : num_max_outEdges = in_degree[i];
	}
	for (int i = 0; i < num_v; i++) {
		if (in_degree[i] == 0) {
			Tree* tree = new Tree();
			TreeNode* currentNode = NULL;
			in_degree_zero_DFS(i, &g, tree, currentNode);
			arr_xiaozhi[num_xiaozhi++] = tree;
		}
	}
	cout << num_xiaozhi << endl;
	//DFS后，确定一个tarjan的开始点u
	int u = 0;
	for(;u< num_v;u++){
		Graph::out_edge_iterator outedgeIt, outedgeEnd;
		tie(outedgeIt, outedgeEnd) = out_edges(u, g);
		Edge e = *outedgeIt;
		property_map<Graph, edge_finished_t>::type edge_finished1 = get(edge_finished, g);
		for (; outedgeIt != outedgeEnd; ++outedgeIt) {
			e = *outedgeIt;
			if (edge_finished1[e])
				continue;
			break;
		}
		if (outedgeIt != outedgeEnd) {
			break;
		}
		else if(edge_finished1[e]){
			break;
		}

	}
	cout << u << endl;
	//Tarjan(&g, u, dfn, low, vis, stack, arr_xiaozhi);

	//加入顶点标签vertex_name_t
	//infile.open("C:\\Users\\xslwyz\\Downloads\\cit-HepTh\\Cit-HepTh-dates.txt", ios::in);

}

void Tarjan(Graph * g,int u, int* dfn, int* low, bool* vis, int* stack,Tree ** arr_xiaozhi) {
	static int dfs_num = 0;
	static int top = 0;
	//DFN[ i ] : 在DFS中该节点被搜索的次序(时间戳)
	dfn[u] = ++dfs_num;
	//LOW[ i ] : 为i或i的子树能够追溯到的最早的栈中节点的次序号
	//当DFN[ i ]==LOW[ i ]时，为i或i的子树可以构成一个强连通分量。
	low[u] = dfs_num;
	vis[u] = true;//是否在栈中
	stack[++top] = u;
	Graph::out_edge_iterator outedgeIt, outedgeEnd;
	tie(outedgeIt, outedgeEnd) = out_edges(0, *g);
	int temp = 0;
	for (; outedgeIt != outedgeEnd; ++outedgeIt)
	{
		Edge e = *outedgeIt;
		Vertex v_tar = target(e, *g);
		temp = v_tar;
		if (!vis[temp]) {
			Tarjan(g, u, dfn, low, vis, stack, arr_xiaozhi);
			low[u] = min(low[u], low[temp]);
		}
		else if (vis[temp])low[u] = min(low[u], dfn[temp]);
	}
	//当dfn[i]==low[i]时，为i或i的子树可以构成一个强连通分量。
	if (dfn[u] == low[u]) {//构成强连通分量
		//vis[u] = false;
		////这里要建小枝模式
		////T={VT,ET,ST,Lable}
		//Tree *tree=new Tree();
		////color[u] = ++col_num;//染色
		//while (stack[top] != u) {//清空
		//	//add_edge(stack[top - 1], stack[top],*tree);
		//	vis[stack[top--]] = false;
		//}
		//if (num_vertices(*gg)) {
		//	arr_xiaozhi[num_xiaozhi] = gg;
		//	num_xiaozhi++;
		//}
		//top--;
	}
}



int main(int argc, char* argv[]){
	preProcess_Cit_Hepth();
	graph_proc();
	return 0;
}

