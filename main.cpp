#include "head.h"

string loc1 = "D:\\cit-HepTh\\Cit-HepTh";
string loc2 = ".txt";
string loc = loc1 + loc2;
string loc_fixed = loc1 + "-fixed" + loc2;
string loc_dic = loc1 + "-dic" + loc2;
string loc_log = loc1 + "-log" + loc2;
static int num_xiaozhi = 0;//小枝的数量
static int num_max_outEdges = 0;
static int num_v = 0;//图的顶点数
typedef property<vertex_name_t, int> VertexProperties;
typedef property<edge_finished_t, int> EdgeProperties;
typedef adjacency_list<listS, vecS, directedS, VertexProperties, EdgeProperties> Graph;
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
	Tree(int val) { root = current = new TreeNode(val); }
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

void Tree::setRoot(int val)
{
	root = current = new TreeNode(val);  //建立根结点
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

//对入度为零的点，进行DFS处理
void in_degree_zero_DFS(int u, Graph& g, Tree* tree, TreeNode* currentNode, vector <int>& in_degree, queue <int>& in_degree_zero_queue) {
	Graph::out_edge_iterator outedgeIt, outedgeEnd,iterTemp;
	tie(outedgeIt, outedgeEnd) = out_edges(u, g);
	Vertex v;
	while (outedgeIt != outedgeEnd) {
		iterTemp = outedgeIt++;
		v = target(*iterTemp, g);
		if (tree->Find(v)) {
			continue;
		}
		TreeNode* treeNode = new TreeNode(v);
		TreeNode* temp = currentNode->left;
		currentNode->left = treeNode;
		treeNode->right = temp;
		remove_edge(*iterTemp, g);
		in_degree[v] = in_degree[v] - 1;//这里要判断删去边后，这个顶点的入度是否变为0，如果变为0，则需要将其加入in_degree_zero_queue
		if (in_degree[v] == 0) {
			in_degree_zero_queue.push(v);
		}
		in_degree_zero_DFS(v, g, tree, currentNode, in_degree, in_degree_zero_queue);
	}
}
void graph_dfs(Graph & g,vector <Tree *> & xiaozhi, vector <int> & in_degree,queue <int> & in_degree_zero_queue) {
	while (in_degree_zero_queue.size()) {//将队列里元素都处理完，就退出并进行tarjan
		int u = in_degree_zero_queue.front();
		if (out_edges(u, g).first == out_edges(u, g).second) {//如果没出边了才POP出，不然会出现提前结束
			in_degree_zero_queue.pop();//取第一个，并pop出
			continue;
		}
		Tree* tree = new Tree(u);//以头元素为根，建树
		TreeNode* currentNode = tree->current;
		in_degree_zero_DFS(u, g, tree, currentNode, in_degree,in_degree_zero_queue);//当前指针指向根，然后进入循环方法in_degree_zero_DFS
		xiaozhi.push_back(tree);//处理完后，将树加入
	}	
}

void graph_tarjan(Graph & g, vector<Tree*>& xiaozhi, int u, std::unordered_map<int, int>& dfn,
	std::unordered_map<int, int>& low,	stack<int>& sta, 
	std::unordered_set<int>& in_stack, queue <int>& in_degree_zero_queue, int dfs_num = 0) {
	Vertex v;
	dfn.insert(pair<int, int>(u, ++dfs_num));//DFN:在DFS中该节点被搜索的次序
	low.insert(pair<int, int>(u, dfs_num));//LOW:为i或i的子树能够追溯到的最早的栈中节点的次序号
	//当DFN[ i ]==LOW[ i ]时，为i或i的子树可以构成一个强连通分量。
	sta.push(u);
	in_stack.insert(u);
	Graph::out_edge_iterator outedgeIt, outedgeEnd,iterTemp;
	tie(outedgeIt, outedgeEnd) = out_edges(u, g);
	while (outedgeIt != outedgeEnd) {
		iterTemp = outedgeIt++;
		v = target(*iterTemp, g);
		if (!dfn.count(v)) {//如果未被访问过就不会在DFN中有值
			graph_tarjan(g, xiaozhi, v, dfn, low, sta, in_stack,in_degree_zero_queue,dfs_num);
			low[u] = min(low[u], low[v]);
		}
		else if (in_stack.count(v)) {//如果v在栈中
			low[u] = min(low[u], dfn[v]);
		}
	}
	if (dfn[u] == low[u]) {
		do {
			v = sta.top();
			sta.pop();
			in_stack.erase(v);
		} while (v != u);
		in_degree_zero_queue.push(u);
	}

}
//G=(V,E,L_V)
//T={V,E,S_T(T中顶点的数量),L_V}
//使用bidirectionalS 替换directedS，这样可以使用in_edge、in_degree等，但消耗要增一倍
void graph_proc() {
	Graph g;
	int x, y = 0;
	string str;
	ifstream infile;
	infile.open(loc_fixed, ios::in);//加入顶点和边
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {//无视文件末尾空行
		istringstream strstream(str);
		strstream >> x >> y;
		if(x!=y)
			add_edge(x, y, g);
	}
	infile.close();   //关闭文件
	num_v = num_vertices(g);

	vector<int> in_degree(num_v, 0);//vector存所有顶点的入度，初始化为0
	infile.open(loc_fixed, ios::in);//再读取一次，因为第一次读取后，才知道顶点个数num_v
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {
		istringstream strstream(str);
		strstream >> x >> y;
		in_degree[y] = in_degree[y] + 1;
	}
	infile.close();   //关闭文件
	queue<int> in_degree_zero_queue;//统计入度是0的所有顶点
	for (int i = 0; i < num_v; i++) {
		if (in_degree[i] == 0)
			in_degree_zero_queue.push(i);
	}
	vector<Tree *> xiaozhi;//存储获得的小枝
	int tarjan_begin_i = 0;//第一个in_degree!=0的点的位置
	cout << num_vertices(g) << endl;
	while (num_edges(g)) {//只要还有边，就重复1.DFS，2.tarjan
		cout << num_edges(g) << endl;
		cout << xiaozhi.size() << endl;
		graph_dfs(g, xiaozhi,in_degree, in_degree_zero_queue);
		//cout << num_edges(g) << endl;
		//cout << xiaozhi.size() << endl;
		std::unordered_map<int, int> dfn;
		std::unordered_map<int, int> low;
		stack<int> sta;
		std::unordered_set<int> in_stack;
		while (tarjan_begin_i < num_v) {//DFS结束了，in_degree=0的点都遍历完了
			if (in_degree[tarjan_begin_i] > 0) {//确定一个tarjan开始的点tarjan_begin_i
				//tarjan TODO
				graph_tarjan(g, xiaozhi, tarjan_begin_i, dfn, low, sta,in_stack, in_degree_zero_queue);
				tarjan_begin_i++;//这个点处理完之后，他的出度一定为0了，没有价值了
				break;//完成后直接退出循环，然后DFS
			}
			tarjan_begin_i++;
		}
		//cout << tarjan_begin_i << endl;
	}
	cout << num_edges(g) << endl;
	cout << xiaozhi.size() << endl;
	ofstream outfile_log;
	outfile_log.open(loc_log, ios::out | ios::trunc);
	Graph::edge_iterator iter1, iter2;
	tie(iter1,iter2)=edges(g);
	while (iter1 != iter2) {
		outfile_log << source(*iter1,g) << '\t' << target(*iter1,g) << endl;
		iter1++;
	}
	outfile_log.close();
}

void stack_tc(Graph& g,int u) {





}
void test() {
	//STACK_TC
	Graph g;
	int x, y = 0;
	string str;
	ifstream infile;
	infile.open(loc_fixed, ios::in);//加入顶点和边
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {//无视文件末尾空行
		istringstream strstream(str);
		strstream >> x >> y;
		if (x != y)
			add_edge(x, y, g);
	}
	infile.close();   //关闭文件
	num_v = num_vertices(g);
	int u = 0;
	vector<int> vis(num_v,0);
	while (u < num_v) {
		if(vis[u]==0)
		stack_tc(g, u);
	}






}

int main(int argc, char* argv[]){
	//preProcess_Cit_Hepth();
	//graph_proc();
	test();
	return 0;
}

