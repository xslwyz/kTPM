#include "head.h"

string loc1 = "C:\\Users\\xslwyz\\Downloads\\cit-HepTh\\Cit-HepTh";
string loc2 = ".txt";
string loc_fixed = loc1 + "-fixed" + loc2;
string loc_dic = loc1 + "-dic" + loc2;


//�����ֵ�
void preProcess_Cit_Hepth() {
	ifstream infile;
	//���붥��ͱ�
	infile.open(loc1+loc2, ios::in);
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

void graph_proc() {
	//G=(V,E,L_V)
	//T={V,E,S_T(T�ж��������),L_V}
	// Property types
	typedef property<vertex_name_t, int, property<vertex_index_t, int, property<vertex_index1_t, int> > > VertexProperties;
	// Graph type
	typedef adjacency_list<vecS, vecS, directedS, VertexProperties> Graph;
	// Graph instance
	Graph g;
	Graph* g_xiaozhi;
	// Property accessors
	property_map<Graph, vertex_name_t>::type vertex_lable = get(vertex_name, g);
	property_map<Graph, vertex_index_t>::type vertex_DFN = get(vertex_index, g);
	property_map<Graph, vertex_index1_t>::type vertex_LOW = get(vertex_index1, g);
	// Create the vertices
	typedef graph_traits<Graph>::edge_descriptor Edge;
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	//Vertex u1;
	//u1 = add_vertex(g);
	//vertex_lable[u1] = 1;
	//Vertex u2;
	//u2 = add_vertex(g);
	//vertex_lable[u2] = 2;
	//Vertex u3;
	//u3 = add_vertex(g);
	//vertex_lable[u3] = 3;
	// Create the edges

	//Edge e1;
	//e1 = (add_edge(u1, u2, g)).first;
	//Edge e2;
	//e2 = add_edge(u1, u3, g).first;

	ifstream infile;
	//���붥��ͱ�
	infile.open(loc_fixed, ios::in);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	int x, y = 0;
	string str;
	while (getline(infile, str)) {            // ��δ���ļ�����һֱѭ��
		istringstream strstream(str);
		strstream >> x >> y;
		add_edge(x, y, g);
	}
	infile.close();   //�ر��ļ�
	int num_v = num_vertices(g);
	//int* DFN;
	////��ʱ�����޳�ֵ
	//DFN = new int[num_v];
	//int* Low;
	//Low = new int[num_v];
	//bool* vis;
	//vis = new bool[num_v];
	//int* stack;
	//stack = new int[num_v];
	Graph::out_edge_iterator outedgeIt, outedgeEnd;
	tie(outedgeIt, outedgeEnd) = out_edges(0, g);
	for (; outedgeIt != outedgeEnd; ++outedgeIt)
	{
		Edge e= *outedgeIt ;
		Vertex v = target(e,g);
	}




	//���붥���ǩvertex_name_t
	//infile.open("C:\\Users\\xslwyz\\Downloads\\cit-HepTh\\Cit-HepTh-dates.txt", ios::in);

}

typedef property<vertex_name_t, int, property<vertex_index_t, int, property<vertex_index1_t, int> > > VertexProperties;
typedef adjacency_list<vecS, vecS, directedS, VertexProperties> Graph;
void Tarjan(Graph * g,int u, int* dfn, int* low, bool* vis, int* stack) {
	static int dfs_num = 0;
	static int top = 0;
	//DFN[ i ] : ��DFS�иýڵ㱻�����Ĵ���(ʱ���)
	dfn[u] = ++dfs_num;
	//LOW[ i ] : Ϊi��i�������ܹ�׷�ݵ��������ջ�нڵ�Ĵ����
	//��DFN[ i ]==LOW[ i ]ʱ��Ϊi��i���������Թ���һ��ǿ��ͨ������
	low[u] = dfs_num;
	vis[u] = true;//�Ƿ���ջ��
	stack[++top] = u;
	typedef graph_traits<Graph>::edge_descriptor Edge;
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	Graph::out_edge_iterator outedgeIt, outedgeEnd;
	tie(outedgeIt, outedgeEnd) = out_edges(0, *g);
	int temp = 0;
	for (; outedgeIt != outedgeEnd; ++outedgeIt)
	{
		Edge e = *outedgeIt;
		Vertex v_tar = target(e, *g);
		temp = v_tar;
		if (!vis[temp]) {
			Tarjan(g, x, dfn, low, vis, stack);
			low[u] = min(low[u], low[temp]);
		}
		else if (vis[temp])low[u] = min(low[u], dfn[temp]);
	}
	//��dfn[i]==low[i]ʱ��Ϊi��i���������Թ���һ��ǿ��ͨ������
	if (dfn[u] == low[u]) {//����ǿ��ͨ����
		vis[u] = false;
		//����Ҫ��С֦ģʽ
		//T={VT,ET,ST,Lable}
		
		//color[u] = ++col_num;//Ⱦɫ
		while (stack[top] != u) {//���
			

			color[stack[top]] = col_num;
			vis[stack[top--]] = false;
		}
		top--;
	}
}



int main(int argc, char* argv[]){
	//preProcess_Cit_Hepth();
	graph_proc();
	return 0;
}

