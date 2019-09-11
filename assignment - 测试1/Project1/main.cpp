/************************************************************************/
/* $Id: MainP.cpp 65 2010-09-08 06:48:36Z yan.qi.asu $                                                                 */
/************************************************************************/

#include <limits>
#include <set>
#include <map>
#include <queue>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include<iterator>
#include "GraphElements.h"
#include "Graph.h"
#include "DijkstraShortestPathAlg.h"
#include "YenTopKShortestPathsAlg.h"
#include<regex>
#include <iomanip>
#include <io.h>

using namespace std;


void testDijkstraGraph()
{
	Graph* my_graph_pt = new Graph("geom_TmCar_combine.txt");
	DijkstraShortestPathAlg shortest_path_alg(my_graph_pt);
	BasePath* result =
		shortest_path_alg.get_shortest_path(
			my_graph_pt->get_vertex(5630), my_graph_pt->get_vertex(6949));
	result->PrintOut(cout);
}

void testYenAlg()
{
	Graph my_graph("geom_Tm_combine.txt");

	YenTopKShortestPathsAlg yenAlg(my_graph, my_graph.get_vertex(1),
		my_graph.get_vertex(33));

	int i = 0;
	while (yenAlg.has_next())
	{
		++i;
		yenAlg.next()->PrintOut(cout);
	}
}

//此函数主要用于生成可被Graph读取的文件，即第一行点个数，剩下行边（起点，终点，权值）的形式
void inputFileConvert(string whichTmDsgnLink) 
{
	const char* file_name = whichTmDsgnLink.c_str();
	string vehi;
	string filename = whichTmDsgnLink;
	string::iterator itB = filename.begin();
	string::iterator itE = filename.end();
	vehi.assign(itB + 10, itE - 4);
	//1. Check the validity of the file
	ifstream ifs_geom("geom.txt");
	ofstream ofs_geom_Tm_combine("geom_Tm"+vehi+"_combine.txt");
	ifstream ifs_Tm(whichTmDsgnLink);
	if (!ifs_geom)
	{
		cerr << "The file " << file_name << " can not be opened!" << endl;
		exit(1);
	}
	int total_vertex_number = 0;
	set<int> vertexes;
	int tempnode_start;
	int tempnode_stop;
	string end_str;
	string end_str_again;
	double middle_double;
	while (ifs_geom>> tempnode_start)
	{
		ifs_geom >> tempnode_stop;
		vertexes.insert(tempnode_start);
		vertexes.insert(tempnode_stop);
		for (int i = 0;i < 6;i++) {
			ifs_geom >> middle_double;
		}
		ifs_geom >> end_str;
		
	}
	ifs_geom.close();
	ifs_geom.open("geom.txt");
	ofs_geom_Tm_combine << vertexes.size()<<endl<<endl;
	getline(ifs_Tm, end_str);
	while (getline(ifs_Tm, end_str)) {
		ifs_geom >> tempnode_start;
		ifs_geom >> tempnode_stop;
		vertexes.insert(tempnode_start);
		vertexes.insert(tempnode_stop);
		for (int i = 0;i < 6;i++) {
			ifs_geom >> middle_double;
		}
		ifs_geom >> end_str_again;
		ofs_geom_Tm_combine <<tempnode_start<<' '<<tempnode_stop<<' '<< end_str<<endl;
	}
	cout << "total number is" << vertexes.size();

	

	
}

//用于OD相加的结果，有起点，终点，OD量
class ODPath {
public:
	int start;
	int stop;
	double ODvalue;

	ODPath(int sta,int sto) {
		start = sta;
		stop = sto;
		ODvalue = 0.0;
	}
	
	void addOD(double OD)
	{
		ODvalue += OD;
	}

	int getStart() {
		return start;
	}

	int getStop() {
		return stop;
	}

	double getOD() {
		return ODvalue;
	}

	bool operator<(const ODPath& p) const //注意这里的两个const
	{
		return (start < p.start);
	}
	bool operator == (const ODPath& x)
	{ 
		return (this->start == x.start) && (this->stop == x.stop); 
	}
};

//用于存储单个中间节点图表（TVolNode），叠加OD，标准输出
class CrossNode {
public:
	int nodeID;
	vector<int> adjacent_nodes;
	int number_of_nodes;
	vector<vector<double>> OD_matrix;
	CrossNode(int nodeID, vector<int> adjacent_nodes) {
		this->nodeID = nodeID;
		this->adjacent_nodes = adjacent_nodes;
		number_of_nodes = adjacent_nodes.size();
		generate_OD_matrix();

	}

	bool operator<(const CrossNode& p) const //注意这里的两个const
	{
		return (nodeID < p.nodeID);
	}
	bool operator == (const CrossNode& x)
	{
		return (this->nodeID == x.nodeID);
	}

	double rowFlow(int row) {
		double result = 0;
		for (int i = 0;i < number_of_nodes;i++) {
			result += OD_matrix[row][i];
		}
		return result;
	}

	double columnFlow(int column) {
		double result = 0;
		for (int i = 0;i < number_of_nodes;i++) {
			result += OD_matrix[i][column];
		}
		return result;
	}

	void add_OD_To_Matrix(int central, int input_node, int output_node, double ODvalue) {
		if (central != nodeID) {
			cerr << "流图中心点与当前结点不符";
		}
		else {
			int input_node_position, output_node_position;
			vector <int>::iterator input_node_it = find(adjacent_nodes.begin(),
				adjacent_nodes.end(), input_node);
			if (input_node_it != adjacent_nodes.end())
			{
				input_node_position = distance(adjacent_nodes.begin(), input_node_it);
			}

			vector <int>::iterator output_node_it = find(adjacent_nodes.begin(),
				adjacent_nodes.end(), output_node);
			if (output_node_it != adjacent_nodes.end())
			{
				output_node_position = distance(adjacent_nodes.begin(), output_node_it);
			}

			OD_matrix[input_node_position][output_node_position] += ODvalue;
		}
	}

	/*void PrintOut(std::ostream& out_stream)
	{
		out_stream << " 交叉口编号： No." << nodeID << endl;
		for (int i = 0; i < (2 + number_of_nodes) * 10 - 1; i++) {
			out_stream << "-";
		}
		out_stream << endl;
		out_stream << setw(9) << "进口-出口";
		for (int i = 0; i < number_of_nodes; i++) {
			out_stream << setw(9) << adjacent_nodes[i];
		}
		out_stream << setw(9) << "进口流量" << endl;
		for (int i = 0; i < (2 + number_of_nodes) * 10 - 1; i++) {
			out_stream << "-";
		}
		out_stream << endl;
		for (int i = 0; i < number_of_nodes; i++) {
			out_stream << setw(9) << adjacent_nodes[i];
			for (int j = 0; j < number_of_nodes; j++) {
				out_stream << setw(9) << OD_matrix[i][j];
			}
			out_stream << setw(9) << rowFlow(i) << endl;
		}
		for (int i = 0; i < (2 + number_of_nodes) * 10 - 1; i++) {
			out_stream << "-";
		}
		out_stream << endl;
		out_stream << setw(9) << "出口流量";
		double total_flow = 0;
		for (int i = 0; i < number_of_nodes; i++) {
			double colume_flow = columnFlow(i);
			out_stream << setw(9) << colume_flow;
			total_flow += colume_flow;
		}
		out_stream << setw(9) << total_flow << endl;
		for (int i = 0; i < (2 + number_of_nodes) * 10 - 1; i++) {
			out_stream << "-";
		}
		out_stream << endl;
		out_stream << endl;

	}*/

	void PrintOut(std::ostream& out_stream)
	{
		out_stream << "\t" << nodeID;
		for (int i = 0;i < number_of_nodes;i++) {
			out_stream << "\t" << adjacent_nodes[i];
		}
		out_stream << endl;
		for (int i = 0;i < number_of_nodes;i++) {
			out_stream << "\t" << adjacent_nodes[i];
			for (int j = 0;j < number_of_nodes;j++) {
				out_stream << "\t" << OD_matrix[i][j];
			}
			out_stream << endl;
		}
	}

private:
	void generate_OD_matrix() {
		int length_of_nodes = adjacent_nodes.size();
		vector<double> matrix_row;
		for (int i = 0;i < length_of_nodes;i++) {
			for (int j = 0;j < length_of_nodes;j++) {
				matrix_row.push_back(0.0);
			}
			OD_matrix.push_back(matrix_row);
			matrix_row.clear();

		}
	}
};


//分配OD，并产生相应文件输出
void assignOD(string whichTmDsgnLink) 
{
	inputFileConvert(whichTmDsgnLink);

	string vehi;
	string filename = whichTmDsgnLink;
	string::iterator itB = filename.begin();
	string::iterator itE = filename.end();
	vehi.assign(itB + 10, itE - 4);

	//加载小区和关键节点的map
	map<int, int> Chnn_map;
	vector<int>critical_key_list;
	ifstream ifs_Chnn("Chnn.txt");
	int Chnn_firstline, Chnn_critical, Chnn_area;
	ifs_Chnn >> Chnn_firstline;
	while (ifs_Chnn >> Chnn_area) {
		ifs_Chnn >> Chnn_critical;
		critical_key_list.push_back(Chnn_critical);
		Chnn_map.insert(make_pair(Chnn_critical, Chnn_area));
	}

	cout << "小区与关键节点加载完毕"<<endl;

	//初始化存储所有加载ODpathes结果的vector
	vector<ODPath>ODpathes;
	ifstream ifs_geom("geom.txt");
	int tempnode_start;
	int tempnode_stop;
	string end_str;
	double middle_double;
	while (ifs_geom >> tempnode_start)
	{
		ifs_geom >> tempnode_stop;
		ODpathes.push_back(ODPath(tempnode_start, tempnode_stop));
		for (int i = 0;i < 6;i++) {
			ifs_geom >> middle_double;
		}
		ifs_geom >> end_str;

	}

	cout << "存储结果的ODPathes生成完毕" << endl;

	//加载OD矩阵
	vector<vector<double>> ODvector;
	vector<double>OD_one_row;
	ifstream ifs_OD("OD"+vehi+".tsd");
	string OD_first_line;
	getline(ifs_OD, OD_first_line);
	double each_OD_value;
	for (int i = 0;i < 23;i++) {
		for (int j = 0;j < 23;j++) {
			ifs_OD >> each_OD_value;
			OD_one_row.push_back(each_OD_value);
		}
		ODvector.push_back(OD_one_row);
		OD_one_row.clear();

	}

	cout << "OD898979矩阵加载完毕" << endl;

	//加载conn.txt数据
	vector<CrossNode> conn_nodes;
	ifstream ifs_conn("conn_convert.txt");
	int conn_id, conn_count, conn_adjacent;
	ifs_conn >> conn_id;//跳过第一行
	while (ifs_conn >> conn_id) {
		ifs_conn >> conn_count;
		vector<int> conn_adjacents;
		for (int i = 0;i < conn_count;i++) {
			ifs_conn >> conn_adjacent;
			conn_adjacents.push_back(conn_adjacent);
		}
		if (conn_count > 2) {
			conn_nodes.push_back(CrossNode(conn_id, conn_adjacents));
		}
	}

	cout << "conn结果矩阵加载完毕" << endl;


	//遍历所有可能的小区组合并计算
	cout << "conn结果矩阵加载完毕" << endl;
	Graph* my_graph_pt = new Graph("geom_Tm"+vehi+"_combine.txt");
	DijkstraShortestPathAlg shortest_path_alg(my_graph_pt);
	cout << "conn结果矩阵加载完毕" << endl;
	for (int i = 0;i < critical_key_list.size();i++) {
		for (int j = i + 1;j < critical_key_list.size();j++) {
			if (i == 0 && j == 18) {
				cout << "asd";
			}
			//计算最短路径
			BasePath* result =
				shortest_path_alg.get_shortest_path(
					my_graph_pt->get_vertex(critical_key_list[i]), my_graph_pt->get_vertex(critical_key_list[j]));

			//根据关键节点获得OD值
			int area1 = Chnn_map.find(critical_key_list[i])->second;
			int area2 = Chnn_map.find(critical_key_list[j])->second;

			/*int size1 = ODvector.size();
			if (area1 - 1 > size1 - 1)
			{
				int size2 = ODvector[area1 - 1].size();
				if (area2 - 1 > size2 - 1)
					int aaa = 5;
			}*/
			
			double current_OD = ODvector[area1 - 1][area2 - 1];

			//将OD值加载到对应路段上
			vector<BaseVertex*> shortest_path_vertex_list = result->getVertexList();
			if (shortest_path_vertex_list.size() != 0) {
				for (int k = 0; k < (shortest_path_vertex_list.size() - 1); k++) {
					vector<ODPath>::iterator rst = find(ODpathes.begin(), ODpathes.end(), ODPath(shortest_path_vertex_list[k]->getID(), shortest_path_vertex_list[k + 1]->getID()));
					rst->addOD(current_OD);
					if (k > 0) {//此处用于生成TVolNode文件，目前因为conn文件尚存在问题无法正常执行。
						vector<CrossNode>::iterator conn_rst = find(conn_nodes.begin(), conn_nodes.end(), CrossNode(shortest_path_vertex_list[k]->getID(), vector<int>()));
						if (conn_rst != conn_nodes.end()) {//中间节点有可能是只连接两个方向的，不一定有对应的表
							//conn_rst->add_OD_To_Matrix(shortest_path_vertex_list[k]->getID(), shortest_path_vertex_list[k - 1]->getID(), shortest_path_vertex_list[k + 1]->getID(), current_OD);
						}

					}
				}
			}


			//重复上述步骤但相反
			//计算最短路径
			BasePath* result_rev =
				shortest_path_alg.get_shortest_path(
					my_graph_pt->get_vertex(critical_key_list[j]), my_graph_pt->get_vertex(critical_key_list[i]));

			//根据关键节点获得OD值
			int area1_rev = Chnn_map.find(critical_key_list[j])->second;
			int area2_rev = Chnn_map.find(critical_key_list[i])->second;
			double current_OD_rev = ODvector[area1_rev - 1][area2_rev - 1];

			//将OD值加载到对应路段上
			vector<BaseVertex*> shortest_path_vertex_list_rev = result_rev->getVertexList();
			if (shortest_path_vertex_list_rev.size() != 0) {
				for (int k = 0; k < shortest_path_vertex_list_rev.size() - 1; k++) {
					vector<ODPath>::iterator rst_rev = find(ODpathes.begin(), ODpathes.end(), ODPath(shortest_path_vertex_list_rev[k]->getID(), shortest_path_vertex_list_rev[k + 1]->getID()));
					rst_rev->addOD(current_OD_rev);
					if (k > 0) {
						vector<CrossNode>::iterator conn_rst = find(conn_nodes.begin(), conn_nodes.end(), CrossNode(shortest_path_vertex_list_rev[k]->getID(), vector<int>()));
						if (conn_rst != conn_nodes.end()) {//中间节点有可能是只连接两个方向的，不一定有对应的表
							//conn_rst->add_OD_To_Matrix(shortest_path_vertex_list_rev[k]->getID(), shortest_path_vertex_list_rev[k - 1]->getID(), shortest_path_vertex_list_rev[k + 1]->getID(), current_OD_rev);
						}

					}
				}
			}
		}
	}
	cout << "OD结果相加完毕,conn结果生成完毕" << endl;

	
	ofstream ofs_geom_Tm_combine("VolLink" + vehi + ".tsd");
	ofs_geom_Tm_combine << "{'format':'table', 'row-header-type':'none', 'column-header-type':'none', 'table-data-type':'float'}" << endl;
	for (int i = 0;i < ODpathes.size();i++) {
		ofs_geom_Tm_combine << ODpathes[i].getOD() << endl;
	}

	ofstream ofs_TVolNode("TVolNode" + vehi + ".tsd");
	//ofs_TVolNode << "{'format':'txt'}" << endl;
	for (int i = 0;i < conn_nodes.size();i++) {
		conn_nodes[i].PrintOut(ofs_TVolNode);
	}

	
}

//获取某一路径下所有文件
void getFiles(string path, vector<string>& files)
{
	//文件句柄  
	long   hFile = 0;
	//文件信息  
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//如果是目录,迭代之  
			//如果不是,加入列表  
			if ((fileinfo.attrib & _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}

class ConnStatistic {
public:
	int central_node_id;
	set<int> adjacents;
	
	ConnStatistic(int central) {
		central_node_id = central;
	}

	void add_adjacent(int adjacent) {
		adjacents.insert(adjacent);
	}

	set<int> getAdjacents() {
		return adjacents;
	}

	bool operator<(const ConnStatistic& p) const //注意这里的两个const
	{
		return (central_node_id < p.central_node_id);
	}
	bool operator == (const ConnStatistic& x)
	{
		return (this->central_node_id == x.central_node_id);
	}
};

void addConnPath(int central, int adjacent,vector<ConnStatistic>&conn_convert ) {
	vector<ConnStatistic>::iterator rst = find(conn_convert.begin(), conn_convert.end(), ConnStatistic(central));
	if (rst != conn_convert.end()) {
		rst->add_adjacent(adjacent);
	}
	else {
		ConnStatistic new_conn_node(central);
		new_conn_node.add_adjacent(adjacent);
		conn_convert.push_back(new_conn_node);
	}
}

//用于将conn文件转化成单行进出节点都包含的形式
void connConvert() {
	vector<ConnStatistic> conn_nodes;
	ifstream ifs_conn("conn.txt");
	int conn_id, conn_count, conn_adjacent;
	ifs_conn >> conn_id;//跳过第一行
	while (ifs_conn >> conn_id) {
		ifs_conn >> conn_count;
		for (int i = 0;i < conn_count;i++) {
			ifs_conn >> conn_adjacent;
			addConnPath(conn_id, conn_adjacent, conn_nodes);
			addConnPath(conn_adjacent,conn_id, conn_nodes);
		}
	}
	sort(conn_nodes.begin(), conn_nodes.end());
	ofstream ofs_conn_convert("conn_convert.txt");
	ofs_conn_convert << conn_nodes.size() << endl;
	for (int i = 0;i < conn_nodes.size();i++) {
		ofs_conn_convert << conn_nodes[i].central_node_id<<" ";
		set<int>adjacents = conn_nodes[i].getAdjacents();
		ofs_conn_convert << adjacents.size()<<" ";
		for (set<int>::iterator it = adjacents.begin();it != adjacents.end();it++)
		{
			ofs_conn_convert << *it <<" ";
		}
		ofs_conn_convert << endl;
	}
}

//筛选出当前路径下所有TmDsgn文件
vector<string> getAllTmDsgnFiles() {
	vector<string>files;
	vector<string>tmdsgnfiles;
	getFiles(".//", files);
	char str[30];
	int size = files.size();
	for (int i = 0;i < size;i++)
	{
		string tmdsgn = "TmDsgn";
		string::size_type idx = files[i].find(tmdsgn);//在a中查找b.
		if (idx != string::npos) {
			tmdsgnfiles.push_back(files[i].substr(4));
		}
	}
	return tmdsgnfiles;
}


int main()
{
	/*cout << "Welcome to the real world!" << endl;
	vector<ODPath>ODpathes;
	ODpathes.push_back(ODPath(1, 3));
	vector<ODPath>::iterator rst = find(ODpathes.begin(), ODpathes.end(), ODPath(1,3));*/
	/*string reguTmstr("//(?<=TmDsgnLink).+?(?=.tsd)");
	regex pattern(reguTmstr, regex::icase);
	smatch result;*/
	/*string vehi;
	string filename = "TmDsgnLinkCar.tsd";
	string::iterator itB = filename.begin();
	string::iterator itE = filename.end();
	vehi.assign(itB+10, itE-4);
	cout << vehi;*/
	/*if (regex_match(reguTmstr, result, pattern)) {
		cout << result[0] << endl;
	}*/
	//testDijkstraGraph();
	//system("pause");
	//inputFileConvert("TmDsgnLinkCar.tsd");
	/*vector<int> uu;
	uu.push_back(5);
	uu.push_back(6);
	cout << uu[0];*/
	//assignOD("TmDsgnLinkCar.tsd");
	//cout.setf(std::ios::right);
	/*cout<< "-------------------------------------------------" << endl;
	cout <<setw(9)<<"进口－出口" <<setw(9)<< 1229<< setw(9) << 33 << setw(9) << 4171 << setw(9) << "进口流量" << endl;
	cout << string("----------------------------------------------------------").size();*/


	//**********************************************重要！！！！！！！！！！！！！*******************************
	/*vector<string>tm_files = getAllTmDsgnFiles(); 这是一个遍历示例，使用时请将需要的Tm文件，OD文件，chnn，conn以及geom等文件全部放在project1文件夹下即可批量生成volLink以及Tvol
	for (int i = 0;i < tm_files.size();i++) {
		assignOD(tm_files[i]);
	}*/
	assignOD("TmDsgnLinkCar.tsd");


	//connConvert();
}
