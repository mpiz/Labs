#include "VFEM_o2_t1_cube.h"

const double eps0 = 8.85418782 * 1e-12;

VFEM_o2_t1_cube:: VFEM_o2_t1_cube() : BaseElement() {
	dofs_n = 1;
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


void VFEM_o2_t1_cube::create_element(vector<int>& nodes_numbs, int el_i) {
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
	const int el_edges_n = 12;
	int el_edges[el_edges_n];
	int el_faces[6];
	vector<dof_type> el_dofs;
	

	static const int edge_nums[el_edges_n][2] = {
		//Вдоль OX
		{0, 1}, {2, 3}, {4, 5}, {6, 7},
		//Вдоль OY
		{0, 2}, {1, 3}, {4, 6}, {5, 7},
		//Вдоль OZ
		{0, 4}, {1, 5}, {2, 6}, {3, 7}
	};

	for(int i = 0; i < el_edges_n; i++)
		el_edges[i] = add_edge(nodes_numbs[edge_nums[i][0]], nodes_numbs[edge_nums[i][1]]);


	// Сформируем грани
	
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
		for(int edge_i = 0; edge_i < el_edges_n; edge_i ++)
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

void VFEM_o2_t1_cube::create_bound_element(vector<int>& nodes_numbs, int el_i) {
	// Упорядочим узлы
	sort(nodes_numbs.begin(), nodes_numbs.end());

	//Знаем, что работает с прямоугольниками - значит всего узлов 4 и считаем, что упорядочены они все как нам надо
	vector<node> el_nodes;
	for(auto iter = nodes_numbs.begin(); iter != nodes_numbs.end(); iter++)
		el_nodes.push_back(nodes[*iter]);

	/*
	==================================
		Нумерация геометрии элемента
	===================================
	*/ 

	// Сформируем все рёбра
	const int el_edges_n = 4;
	int el_edges[el_edges_n];
	int el_face_n;
	vector<dof_type> el_dofs;

	static const int edge_nums[el_edges_n][2] = {
		//Вдоль OY
		{0, 2}, {1, 3},
		//Вдоль OX	
		{0, 1}, {2, 3}
	};

	for(int i = 0; i < el_edges_n; i++)
		el_edges[i] = add_edge(nodes_numbs[edge_nums[i][0]], nodes_numbs[edge_nums[i][1]]);


	// Получим номер грани
	el_face_n = add_face(nodes_numbs[0], nodes_numbs[1], nodes_numbs[2], nodes_numbs[3]);



	//А вот теперь получим степени свободы... ох ёжкин...
	// Сначала 12 базисных функций первого порядка првого типа, затем столько же первого порядка второго типа
	for (int type = 1; type < 3; type++)
		for(int edge_i = 0; edge_i < el_edges_n; edge_i ++)
			el_dofs.push_back(get_geom_dof(GEOM_TYPES::EDGE, el_edges[edge_i], type));

	// Теперь получим 2 функции от грани
		for(int type = 1; type < 3; type++)
			el_dofs.push_back(get_geom_dof(GEOM_TYPES::FACE, el_face_n, type));


	bound_faces[el_i] = quadelement(el_nodes, el_dofs);

	//Запишим все граничные степени свободы
	for(auto dof_i = el_dofs.begin(); dof_i != el_dofs.end(); dof_i++)
		first_bound_dofs_set.insert(*dof_i);

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

void VFEM_o2_t1_cube::input_bound(string file_elements) {

	ifstream inp_bound(file_elements.c_str());

	inp_bound >> first_bound_elements_n;
	bound_faces.resize(first_bound_elements_n);

	vector<int> el_nodes;
	el_nodes.resize(4);
	for(int el_i = 0; el_i < first_bound_elements_n; el_i++) {
		for(int i = 0; i < 4; i++)
			inp_bound >> el_nodes[i];
		create_bound_element(el_nodes, el_i);

	}

	inp_bound.close();
	
}

void VFEM_o2_t1_cube::prepare_boundry() {
	
	//Перепишем из set'а в вектор
	for(auto dof_i = first_bound_dofs_set.begin(); dof_i != first_bound_dofs_set.end(); dof_i++)
		first_bound_dofs.push_back(*dof_i);
	first_bound_n = first_bound_dofs_set.size();

	boundory_method.set_elements(bound_faces, first_bound_dofs_set);
	first_bound_dofs_set.clear();
	bound_faces.clear();

	boundory_method.calculate();
}

void VFEM_o2_t1_cube::calculate() {
	prepare_boundry();

	vector<vfunc3d> rps;
	rps.push_back(
		[](double x, double y, double z)->vec3d{
			return 1e8*eps0 * vec3d(x, 0, 0);
		}
	);

	generate_port();
	generate_matrix_with_out_bound(rps);

	// Учёт первых краевых
	int basis_i = 0;


#ifdef DEBUGOUTP
	auto rp_v = rp[basis_i];
	auto sol = solutions[basis_i];

	output_matrix("matrix_before.txt");
#endif

	for(int k = 0; k < first_bound_n; k++)	{
		int cur_row = first_bound_dofs[k];
		double val = boundory_method.get_dof_weight(cur_row);

		rp[basis_i][cur_row] = val;

		di[cur_row] = 1;

		int i_s = gi[cur_row], i_e = gi[cur_row+1];
		for(int i = i_s; i < i_e; i++){
			rp[basis_i][gj[i]] -= gg[i]*val;
			gg[i] = 0;
		}
		for(int p = cur_row + 1; p < local_dof_n; p++){
			int i_s = gi[p], i_e = gi[p+1];
			for(int i = i_s; i < i_e; i++){
				if(gj[i] == cur_row){
					rp[basis_i][p] -= gg[i]*val;
					gg[i] = 0;
				}
			}
		}
	}

#ifdef DEBUGOUTP
	output_matrix("matrix_after.txt");
#endif

	solver.init(&gi.front(), &gj.front(), &di.front(), &gg.front(), local_dof_n);
	solver.solve(rp[basis_i], solutions[basis_i]);

	point p1(0.8, 0.8, 0.8);
	auto res = value_element_vec(find_element(p1), p1);
}


void VFEM_o2_t1_cube::output_weights(string file_name, string file_bound) {
	ofstream outp_file(file_name.c_str());

	for(int i = 0; i < local_dof_n; i++)
		outp_file << solutions[0][i] << endl;

	outp_file.close();

	boundory_method.output_weights(file_bound);

}


cubeelement* VFEM_o2_t1_cube::find_element(point pn) {
	for(int i = 0; i < elements_n; i++)
		if (elements[i].in_element(pn)) 
			return &elements[i];

	return nullptr;
}
	
double VFEM_o2_t1_cube::get_lambda(cubeelement& el){
	return 1.0;
}
