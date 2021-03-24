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
class TreeNode  //��ʾ���Ľ�㣬�������ֵ�
{
	friend Tree;
public:
	int data;//����
	TreeNode* left;
	TreeNode* right;  //��Ů���ֵ�ָ��
	TreeNode(int value = 0, TreeNode* l = NULL, TreeNode* r = NULL) :data(value), left(l), right(r) {}//���캯��
};
class Tree
{
public:
	Tree() { root = current = NULL; }; //������
	void setRoot(int value);
	TreeNode* getRoot();
	int toLeft();
	int toRight();
	int Parent();
	int Find(int target);
	TreeNode* root, * current;    //��ʾ���͵�ǰָ��
//void PreOrder(ostream& out, TreeNode* p)
	int Find(TreeNode* p, int target);
	//void RemovesubTree(TreeNode* p)
	int FindParent(TreeNode* t, TreeNode* p);
};

void Tree::setRoot(int value)
{
	root = current = new TreeNode(value);  //���������
}

TreeNode* Tree::getRoot()  //�������Ҹ�,ʹ֮��Ϊ��ǰ���
{
	if (root == NULL) { current = NULL; return NULL; }
	else { current = root; return root; }
}

int Tree::toLeft()
{ //Ѱ�ҵ�ǰ���ĵ�һ����Ů.������Ů,��������0,��ǰָ��ΪNULL;
  //���򷵻�1,��ǰָ���Ƶ���ǰ���ĵ�һ����Ů;
	if (current != NULL && current->left != NULL)
	{
		current = current->left; return 1;
	}
	current = NULL; return 0;
}


int Tree::toRight()
{
	//Ѱ�ҵ�ǰ������һ���ֵ�.������һ���ֵ�,��������0,��ǰָ��ΪNULL;���򷵻�1,��ǰָ���ƶ�����ǰ������һ���ֵ�.
	if (current != NULL && current->right != NULL)
	{
		current = current->right; return 1;
	}
	current = NULL; return 0;
}

int Tree::Parent()
{//Ѱ�ҵ�ǰ����˫�׽��.����˫�׽��,��������0,��ǰָ��ΪNULL;���򷵻�1,��ǰָ���Ƶ���ǰ����˫�׽��
	TreeNode* p = current;//p���浱ǰָ��
	TreeNode* t;
	if (current == NULL || current == root) { current = NULL; return 0; }//û��˫��
	t = root;   //t�Ӹ���㿪ʼ
	int  k = FindParent(t, p);  //�Ӹ����t��������ҽ��p��˫�׽��
	return k;
}

//t����������ָ��root
int Tree::FindParent(TreeNode* t, TreeNode* q)
{//˽�к���:�Ӹ�ָ��t��ָ���������,�ҽ��p��˫�׽��,���� current�� 
	TreeNode* p = t->left;//�ҵ�һ������
	while (q != NULL && q != p)       //q==NULL,������;q==p,�ҵ�˫��
	{
		int i = FindParent(q, p);
		if (i != 0) return i;
		q = q->right; //����һ������
	}
	if (q != NULL && q == p) { current = t; return 1; } //�ɹ�����1
	else return 0;    //���򷵻�0
}

int Tree::Find(int target)
{//�������������ϸ���ֵtarget�Ľ��.�����ɹ�,��������1,���򷵻�0
	if (root == NULL) return 0;      //Ϊ����
	return Find(root, target);    //����˽�к���Find����
}

int Tree::Find(TreeNode* p, int target)
{ //˽�к���:���ɸ�ָ��p��ָ���������������������ݳ�Ա����target�Ľ��.�����ɹ�,��������1,ͬʱ�ý���Ϊ��ǰ���;����������0
	int result = 0;
	if (p->data == target) { result = 1; current = p; } //�����ɹ�,p��Ϊ��ǰ���
	else           //����,��������
	{
		TreeNode* q = p->left;
		while (q != NULL && !(result = Find(q, target))) q = q->right;
	}
	return result;
}




//�����ֵ�
void preProcess_Cit_Hepth() {
	ifstream infile;
	//���붥��ͱ�
	infile.open(loc, ios::in);
	ofstream outfile_fixed;
	ofstream outfile_dic;
	outfile_fixed.open(loc_fixed,ios::out|ios::trunc);
	outfile_dic.open(loc_dic, ios::out | ios::trunc);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	string u, v = "";
	//�ӵ����п�ʼȡ��
	int begin = 5;
	string str;
	map<string, int> dic;
	int index_dic = 0;
	map<string, int>::iterator finder;

	while (getline(infile,str)) {            // ��δ���ļ�����һֱѭ��
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
	infile.close();   //�ر��ļ�
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
		//��currentNode->left���ֵܣ�currentNode->left�ǲ���NULL������ν
		TreeNode* temp = currentNode->left;
		currentNode->left = treeNode;
		treeNode->right = temp;
		//g�б��������
		edge_finished1[e] = 1;
		in_degree_zero_DFS(v_tar, g, tree, currentNode);
	}

}

void graph_proc() {
	
	//G=(V,E,L_V)
	//T={V,E,S_T(T�ж��������),L_V}
	//ʹ��bidirectionalS �滻directedS����������ʹ��in_edge��in_degree�ȣ�������Ҫ��һ��
	Graph g;
	int x, y = 0;
	string str;
	ifstream infile;
	//���붥��ͱ�
	infile.open(loc_fixed, ios::in);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {            // ��δ���ļ�����һֱѭ��
		istringstream strstream(str);
		strstream >> x >> y;
		add_edge(x, y, g);
	}
	infile.close();   //�ر��ļ�
	int num_v = num_vertices(g);
	int* dfn= new int[num_v];
	int* low= new int[num_v];
	bool* vis= new bool[num_v];
	int* stack= new int[num_v];
	Tree** arr_xiaozhi =new Tree *[num_v];

	int* in_degree = new int[num_v];
	for (int i = 0; i < num_v; i++)
		in_degree[i] = 0;
	//���붥��ͱ�
	infile.open(loc_fixed, ios::in);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {            // ��δ���ļ�����һֱѭ��
		istringstream strstream(str);
		strstream >> x >> y;
		in_degree[y] = in_degree[y] + 1;
	}

	infile.close();   //�ر��ļ�
	//Graph::out_edge_iterator outedgeIt, outedgeEnd;
	//tie(outedgeIt, outedgeEnd) = out_edges(0, g);
	//for (; outedgeIt != outedgeEnd; ++outedgeIt)
	//{
	//	Edge e = *outedgeIt;
	//	Vertex v = target(e, g);
	//}

	//�������ж��㣬�����Ϊ0�����ȫ�����
	//��bgl̫���ˣ�ֱ�ӽ������
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
	//DFS��ȷ��һ��tarjan�Ŀ�ʼ��u
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

	//���붥���ǩvertex_name_t
	//infile.open("C:\\Users\\xslwyz\\Downloads\\cit-HepTh\\Cit-HepTh-dates.txt", ios::in);

}

void Tarjan(Graph * g,int u, int* dfn, int* low, bool* vis, int* stack,Tree ** arr_xiaozhi) {
	static int dfs_num = 0;
	static int top = 0;
	//DFN[ i ] : ��DFS�иýڵ㱻�����Ĵ���(ʱ���)
	dfn[u] = ++dfs_num;
	//LOW[ i ] : Ϊi��i�������ܹ�׷�ݵ��������ջ�нڵ�Ĵ����
	//��DFN[ i ]==LOW[ i ]ʱ��Ϊi��i���������Թ���һ��ǿ��ͨ������
	low[u] = dfs_num;
	vis[u] = true;//�Ƿ���ջ��
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
	//��dfn[i]==low[i]ʱ��Ϊi��i���������Թ���һ��ǿ��ͨ������
	if (dfn[u] == low[u]) {//����ǿ��ͨ����
		//vis[u] = false;
		////����Ҫ��С֦ģʽ
		////T={VT,ET,ST,Lable}
		//Tree *tree=new Tree();
		////color[u] = ++col_num;//Ⱦɫ
		//while (stack[top] != u) {//���
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

