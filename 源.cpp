#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <queue>
#include <unordered_set>
#include <stack>
using namespace std;

static int num_v_datagraph = 27770;//图的顶点数
static int num_labletype = 10;//标签类型数
const static int default_weight = 999;
const static int topK = 3;
static int num_twig = 0;//小枝的数量
static int num_max_outEdges = 0;

string loc1 = "D:\\cit-HepTh\\Cit-HepTh";
string loc2 = ".txt";
string loc = loc1 + loc2;
string loc_fixed = loc1 + "-fixed" + loc2;
string loc_dic = loc1 + "-dic" + loc2;
string loc_log = loc1 + "-log" + loc2;


class Vertex {
public:
	int value;
	int lable;
	vector<Vertex*> children;
	vector<vector<Vertex*>> children_withlable;//children[a]带有a标签的顶点的集合
	Vertex(int value) :value(value), lable(value % 10) {}
	Vertex(int value, int lable) :value(value), lable(lable) {}
};
class Edge {
public:
	Vertex* source;
	Vertex* target;
	int weight;
	int del_falg = 0;
	Edge(Vertex* source, Vertex* target) :source(source), target(target), weight(default_weight) {}
	Edge(Vertex* source, Vertex* target, int weight) :source(source), target(target), weight(weight) {}
};

class Graph {
public:
	vector<Vertex*> vertex_vec;
	vector<vector<Edge*>> in_edges;
	vector<vector<Edge*>> out_edges;
	vector<vector<Vertex*>> vertex_vec_withlable;//save vertex* by lable

	Graph():vertex_vec(vector<Vertex*>(num_v_datagraph)), in_edges(vector<vector<Edge*>>(num_v_datagraph)),
		out_edges(vector<vector<Edge*>>(num_v_datagraph)), vertex_vec_withlable(vector<vector<Vertex*>>(num_labletype)) {
		vertex_vec;
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
				add_edge(x, y, default_weight);
		}
		infile.close();   //关闭文件
		num_v_datagraph = vertex_vec.size();
	};
	Vertex* getVertex(int x) {
		return vertex_vec[x];
	}
	int in_degree(int x) {
		return in_edges[x].size();
	}
	int out_degree(int x) {
		return out_edges[x].size();
	}
	void add_edge(int x, int y, int weight);
	//基于顶点标签的预处理
	void vertex_children_withlable_preprocessing();
};
void Graph::add_edge(int x, int y, int weight= default_weight) {
	if (vertex_vec[x] == NULL)
		vertex_vec[x] = new Vertex(x);
	if (vertex_vec[y] == NULL)
		vertex_vec[y] = new Vertex(y);
	Edge* e = new Edge(vertex_vec[x], vertex_vec[y], weight);
	in_edges[y].push_back(e);
	out_edges[x].push_back(e);
}
void Graph::vertex_children_withlable_preprocessing() {
	for (int i = 0; i < num_v_datagraph; i++) {
		auto edges = out_edges[i];
		auto source_vertex = getVertex(i);
		Vertex* target_vertex;
		for (int j = 0; j < edges.size(); j++) {
			target_vertex = edges[j]->target;
			auto tar_lable = target_vertex->lable;
			source_vertex->children_withlable[tar_lable].push_back(target_vertex);
			vertex_vec_withlable[tar_lable].push_back(target_vertex);
		}
	}
}

struct Match {
	vector<Vertex*> vertex_vec;
	int score;
};

vector<Match> matches;

void kTPM() {


}

int main()
{
	Graph* g = new Graph();



	std::cout << "Hello World!\n";
}
