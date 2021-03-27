#include "head.h"

string loc1 = "D:\\cit-HepTh\\Cit-HepTh";
string loc2 = ".txt";
string loc = loc1 + loc2;
string loc_fixed = loc1 + "-fixed" + loc2;
string loc_dic = loc1 + "-dic" + loc2;
string loc_log = loc1 + "-log" + loc2;
static int num_xiaozhi = 0;//С֦������
static int num_max_outEdges = 0;
static int num_v = 0;//ͼ�Ķ�����
typedef property<vertex_name_t, int> VertexProperties;
typedef property<edge_finished_t, int> EdgeProperties;
typedef adjacency_list<listS, vecS, directedS, VertexProperties, EdgeProperties> Graph;
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
	Tree(int val) { root = current = new TreeNode(val); }
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

void Tree::setRoot(int val)
{
	root = current = new TreeNode(val);  //���������
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

//�����Ϊ��ĵ㣬����DFS����
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
		in_degree[v] = in_degree[v] - 1;//����Ҫ�ж�ɾȥ�ߺ�������������Ƿ��Ϊ0�������Ϊ0������Ҫ�������in_degree_zero_queue
		if (in_degree[v] == 0) {
			in_degree_zero_queue.push(v);
		}
		in_degree_zero_DFS(v, g, tree, currentNode, in_degree, in_degree_zero_queue);
	}
}
void graph_dfs(Graph & g,vector <Tree *> & xiaozhi, vector <int> & in_degree,queue <int> & in_degree_zero_queue) {
	while (in_degree_zero_queue.size()) {//��������Ԫ�ض������꣬���˳�������tarjan
		int u = in_degree_zero_queue.front();
		if (out_edges(u, g).first == out_edges(u, g).second) {//���û�����˲�POP������Ȼ�������ǰ����
			in_degree_zero_queue.pop();//ȡ��һ������pop��
			continue;
		}
		Tree* tree = new Tree(u);//��ͷԪ��Ϊ��������
		TreeNode* currentNode = tree->current;
		in_degree_zero_DFS(u, g, tree, currentNode, in_degree,in_degree_zero_queue);//��ǰָ��ָ�����Ȼ�����ѭ������in_degree_zero_DFS
		xiaozhi.push_back(tree);//������󣬽�������
	}	
}

void graph_tarjan(Graph & g, vector<Tree*>& xiaozhi, int u, std::unordered_map<int, int>& dfn,
	std::unordered_map<int, int>& low,	stack<int>& sta, 
	std::unordered_set<int>& in_stack, queue <int>& in_degree_zero_queue, int dfs_num = 0) {
	Vertex v;
	dfn.insert(pair<int, int>(u, ++dfs_num));//DFN:��DFS�иýڵ㱻�����Ĵ���
	low.insert(pair<int, int>(u, dfs_num));//LOW:Ϊi��i�������ܹ�׷�ݵ��������ջ�нڵ�Ĵ����
	//��DFN[ i ]==LOW[ i ]ʱ��Ϊi��i���������Թ���һ��ǿ��ͨ������
	sta.push(u);
	in_stack.insert(u);
	Graph::out_edge_iterator outedgeIt, outedgeEnd,iterTemp;
	tie(outedgeIt, outedgeEnd) = out_edges(u, g);
	while (outedgeIt != outedgeEnd) {
		iterTemp = outedgeIt++;
		v = target(*iterTemp, g);
		if (!dfn.count(v)) {//���δ�����ʹ��Ͳ�����DFN����ֵ
			graph_tarjan(g, xiaozhi, v, dfn, low, sta, in_stack,in_degree_zero_queue,dfs_num);
			low[u] = min(low[u], low[v]);
		}
		else if (in_stack.count(v)) {//���v��ջ��
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
//T={V,E,S_T(T�ж��������),L_V}
//ʹ��bidirectionalS �滻directedS����������ʹ��in_edge��in_degree�ȣ�������Ҫ��һ��
void graph_proc() {
	Graph g;
	int x, y = 0;
	string str;
	ifstream infile;
	infile.open(loc_fixed, ios::in);//���붥��ͱ�
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {//�����ļ�ĩβ����
		istringstream strstream(str);
		strstream >> x >> y;
		if(x!=y)
			add_edge(x, y, g);
	}
	infile.close();   //�ر��ļ�
	num_v = num_vertices(g);

	vector<int> in_degree(num_v, 0);//vector�����ж������ȣ���ʼ��Ϊ0
	infile.open(loc_fixed, ios::in);//�ٶ�ȡһ�Σ���Ϊ��һ�ζ�ȡ�󣬲�֪���������num_v
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {
		istringstream strstream(str);
		strstream >> x >> y;
		in_degree[y] = in_degree[y] + 1;
	}
	infile.close();   //�ر��ļ�
	queue<int> in_degree_zero_queue;//ͳ�������0�����ж���
	for (int i = 0; i < num_v; i++) {
		if (in_degree[i] == 0)
			in_degree_zero_queue.push(i);
	}
	vector<Tree *> xiaozhi;//�洢��õ�С֦
	int tarjan_begin_i = 0;//��һ��in_degree!=0�ĵ��λ��
	cout << num_vertices(g) << endl;
	while (num_edges(g)) {//ֻҪ���бߣ����ظ�1.DFS��2.tarjan
		cout << num_edges(g) << endl;
		cout << xiaozhi.size() << endl;
		graph_dfs(g, xiaozhi,in_degree, in_degree_zero_queue);
		//cout << num_edges(g) << endl;
		//cout << xiaozhi.size() << endl;
		std::unordered_map<int, int> dfn;
		std::unordered_map<int, int> low;
		stack<int> sta;
		std::unordered_set<int> in_stack;
		while (tarjan_begin_i < num_v) {//DFS�����ˣ�in_degree=0�ĵ㶼��������
			if (in_degree[tarjan_begin_i] > 0) {//ȷ��һ��tarjan��ʼ�ĵ�tarjan_begin_i
				//tarjan TODO
				graph_tarjan(g, xiaozhi, tarjan_begin_i, dfn, low, sta,in_stack, in_degree_zero_queue);
				tarjan_begin_i++;//����㴦����֮�����ĳ���һ��Ϊ0�ˣ�û�м�ֵ��
				break;//��ɺ�ֱ���˳�ѭ����Ȼ��DFS
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
	infile.open(loc_fixed, ios::in);//���붥��ͱ�
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	while (getline(infile, str)) {//�����ļ�ĩβ����
		istringstream strstream(str);
		strstream >> x >> y;
		if (x != y)
			add_edge(x, y, g);
	}
	infile.close();   //�ر��ļ�
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

