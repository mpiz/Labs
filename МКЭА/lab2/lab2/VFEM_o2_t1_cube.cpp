#include "VFEM_o2_t1_cube.h"

VFEM_o2_t1_cube:: VFEM_o2_t1_cube() : BaseElement() {
	dofs_n = 0;
	local_dof_n = 0;
	cur_edge_count = 0;
	cur_face_count = 0;
}



dof_type VFEM_o2_t1_cube::get_geom_dof(GEOM_TYPES geom, int n_g, int num) {
	stringstream strs;
	switch(geom) {
		case GEOM_TYPES::EDGE :
			strs << "edge_" ;
			break;
		case GEOM_TYPES::FACE :
			strs << "face_" ;
			break;
		case GEOM_TYPES::VOLUME :
			strs << "volume_" ;
			break;
	}
	strs << n_g << "#" << num;

	string key = strs.str();

	if (dofs_map.find(key) != dofs_map.end())
		return dofs_map[key];
	
	int cdc = local_dof_n;
	dofs_map[key] = cdc;
	local_dof_n++;
	return cdc;
}

int VFEM_o2_t1_cube::add_edge(int n1, int n2) {
	auto key = make_tuple(n1, n2);
	auto iter = edges.find(key);
	if (iter != edges.end())
		return iter->second;

	int cec = cur_edge_count;
	edges[key] = cec;
	cur_edge_count++;
	return cec;
}

int VFEM_o2_t1_cube::add_face(int n1, int n2, int n3, int n4) {
	auto key = make_tuple(n1, n2, n3, n4);
	auto iter = faces.find(key);
	if (iter != faces.end())
		return iter->second;

	int cfc = cur_face_count;
	faces[key] = cfc;
	cur_face_count++;
	return cfc;

}


void VFEM_o2_t1_cube::create_element(vector<int> nodes_numbs, int el_i) {
	// Упорядочим узлы
	sort(nodes_numbs.begin(), nodes_numbs.end());

	//Знаем, что работает с парллелипипедами - значит всего узлов 8 и считаем, что упорядочены они все как нам надо
	vector<node> el_nodes;
	for(auto iter = nodes_numbs.begin(); iter != nodes_numbs.end(); iter++)
		el_nodes.push_back(nodes[*iter]);

	/*
	==================================
		Нумерация геометрии элемента
	===================================
	*/ 

	// Сформируем все рёбра
	int el_edges[12];
	int el_faces[6];
	vector<dof_type> el_dofs;

	static const int edge_nums[12][2] = {
		//Вдоль OX
		{0, 1}, {2, 3}, {4, 5}, {6, 7},
		//Вдоль OY
		{0, 2}, {1, 3}, {4, 6}, {5, 7},
		//Вдоль OZ
		{0, 4}, {1, 5}, {2, 6}, {3, 7}
	};

	for(int i = 0; i < 12; i++)
		el_edges[i] = add_edge(nodes_numbs[edge_nums[i][0]], nodes_numbs[edge_nums[i][1]]);


	// Сформируем грани грани
	
	static const int face_nums[6][4] = {
		// Вдоль YZ
		{0, 2, 4, 6}, {1, 3, 5, 7},
		//Вдоль XZ
		{0, 1, 4, 5}, {2, 3, 6, 7},
		//Вдоль XY
		{0, 1, 2, 3}, {4, 5, 6, 7}
	};

	for(int i = 0; i < 6; i++)
		el_faces[i] = add_face(nodes_numbs[face_nums[i][0]], nodes_numbs[face_nums[i][1]], nodes_numbs[face_nums[i][2]], nodes_numbs[face_nums[i][3]]);

	//А вот теперь получим степени свободы... ох ёжкин...
	// Сначала 12 базисных функций первого порядка првого типа, затем столько же первого порядка второго типа
	for (int type = 1; type < 3; type++)
		for(int edge_i = 0; edge_i < 12; edge_i ++)
			el_dofs.push_back(get_geom_dof(GEOM_TYPES::EDGE, el_edges[edge_i], type));

	// Теперь получим по 2 функции от каждой грани
	for(int face_i = 0; face_i < 6; face_i++)
		for(int type = 1; type < 3; type++)
			el_dofs.push_back(get_geom_dof(GEOM_TYPES::FACE, el_faces[face_i], type));

	//И наконец 3 функции, связанные с самим элементом:
	for(int type = 1; type < 4; type++)
		el_dofs.push_back(get_geom_dof(GEOM_TYPES::VOLUME, el_i, type));

	// Заканчиваем формировать элемент
	elements[el_i] = cubeelement(el_nodes, el_dofs);

}


void VFEM_o2_t1_cube::input_mesh(string file_nodes, string file_elements) {
	ifstream inp_node(file_nodes.c_str());

	inp_node >> nodes_n;
	nodes.resize(nodes_n);
	double x, y, z;
	for(int i = 0; i < nodes_n; i++) {
		inp_node >> x >> y >> z;
		nodes[i] = node(x, y, z);
		nodes[i].number = i;
	}

	inp_node.close();

	ifstream inp_elem(file_elements.c_str());
	vector<int> el_nodes;
	el_nodes.resize(8);

	inp_elem >> elements_n;
	elements.resize(elements_n);

	for(int el_i = 0; el_i < elements_n; el_i++) {
		for(int n_i = 0; n_i < 8; n_i++)
			inp_elem >> el_nodes[n_i];
		create_element(el_nodes, el_i);
	}

	inp_elem.close();


}

double VFEM_o2_t1_cube::get_lambda(cubeelement& el){
	return 1.0;
}