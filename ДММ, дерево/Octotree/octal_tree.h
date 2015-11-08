#pragma once

#include <vector>
#include <array>
#include <cstdio>
#include <cmath>

using namespace std;

/*
Небольшая документация классов. Два шаблона образуют восьмиричное дерево.
Параметры шаблонов и требования к ним:
element_type - тип элемента, должен содежрать две функции:
	bool valid_for_tree_node(double x0, double x1, double y0, double y1, double z0, double z1) - вроеряет надо ли добавлять элемент в парралелипипед
	bool in_element(double x, double y, double z) - проверяет принадлежит ли точка (x,y,z) элементу

*/


template<class element_type> class octal_tree_node {
 public:
	 octal_tree_node();
	 void init_node(double x0_s, double x1_s, double y0_s, double y1_s, double z0_s, double z1_s, int split_constant_s, int level_barier_s);
	 void add_element(vector<element_type*>& added_elements);
	 bool add_element(element_type* added_element);
	 element_type* find_point(double x, double y, double z);
	 int get_level();
	 void set_level(int level_s);
	

 private:

	bool point_in_node(double x, double y, double z);
	bool split_condition();
	void split();

	vector<element_type*> elements;
	bool is_node_terminal;
	octal_tree_node<element_type> *sub_nodes;
	double x0, x1, y0, y1, z0, z1;
	int elements_n;
	int split_constant;
	int level;
	int level_barier;


};

template<class element_type> class octal_tree {
 public:
	 void construct_tree(double x0, double x1, double y0, double y1, double z0, double z1, vector<element_type>& elements); //габариты области и список элементов
	 element_type* find_point(double x, double y, double z);

 private:
	 octal_tree_node<element_type> root_node;
	 int total_n;
	 int split_constant;
};


// ==================================== Реализации ====================================

template<class element_type> octal_tree_node<element_type>::octal_tree_node() {
	is_node_terminal = true;
	elements_n = 0;
}

template<class element_type> int octal_tree_node<element_type>::get_level() {
	return level;
}

template<class element_type> void octal_tree_node<element_type>::set_level(int level_s) {
	level = level_s;
}

template<class element_type> void octal_tree_node<element_type>::init_node(double x0_s, double x1_s, double y0_s, double y1_s, double z0_s, double z1_s, int split_constant_s, int level_barier_s) {
	x0 = x0_s;
	x1 = x1_s;
	y0 = y0_s;
	y1 = y1_s;
	z0 = z0_s;
	z1 = z1_s;
	split_constant = split_constant_s;
	level_barier = level_barier_s;
}

template<class element_type> bool octal_tree_node<element_type>::point_in_node(double x, double y, double z) {
	if(x >= x0 && x <= x1 && y >= y0 && y <= y1 && z >= z0 && z <= z1)
		return true;
	return false;

}

template<class element_type> bool octal_tree_node<element_type>::split_condition() {
	if(elements_n + 1 == split_constant &&  fabs(x1 - x0) > 1e-8 && fabs(y1 - y0) > 1e-8 && fabs(z1 - z0) > 1e-8 && level > level_barier)
		return true;
	return false;

}

template<class element_type> bool octal_tree_node<element_type>::add_element(element_type* added_element) {

	if(added_element->valid_for_tree_node(x0, x1, y0, y1, z0, z1)) {

		if(is_node_terminal) { //если узел терминальный (последний)

			if(split_condition()) { //если узел после добавляения перестанет быть терминальным
				elements.push_back(added_element);
				split();
			}

			else { //если после добавления останется треминальным
				elements.push_back(added_element);
			}
			elements_n++;

		}
		else { //если узел не терминальный (есть поддеревья)
			elements_n++;
			for(int i = 0; i < 8; i++) {
					sub_nodes[i].add_element(added_element);
			}

		}
		return true;
	}
	return false;
}

template<class element_type> void octal_tree_node<element_type>::add_element(vector<element_type*>& added_elements) {

	for(int i = 0 ; i < added_elements.size(); i++) {
		if(i == 42)
			i = i;
		add_element(added_elements[i]);
	}

}

template<class element_type> void octal_tree_node<element_type>::split() {
	is_node_terminal = false;
	sub_nodes = new octal_tree_node<element_type> [8];

	double dx = (x1 - x0) / 2.0, dy = (y1 - y0) / 2.0, dz = (z1 - z0) / 2.0;

	sub_nodes[0].init_node(x0, x0+dx, y0, y0+dy, z0, z0+dz, split_constant, level_barier);
	sub_nodes[1].init_node(x0+dx, x1, y0, y0+dy, z0, z0+dz, split_constant, level_barier);
	sub_nodes[2].init_node(x0, x0+dx, y0+dy, y1, z0, z0+dz, split_constant, level_barier);
	sub_nodes[3].init_node(x0+dx, x1, y0+dy, y1, z0, z0+dz, split_constant, level_barier);

	sub_nodes[4].init_node(x0, x0+dx, y0, y0+dy, z0+dz, z1, split_constant, level_barier);
	sub_nodes[5].init_node(x0+dx, x1, y0, y0+dy, z0+dz, z1, split_constant, level_barier);
	sub_nodes[6].init_node(x0, x0+dx, y0+dy, y1, z0+dz, z1, split_constant, level_barier);
	sub_nodes[7].init_node(x0+dx, x1, y0+dy, y1, z0+dz, z1, split_constant, level_barier);

	for(int i = 0; i < 8; i++)
		sub_nodes[i].set_level(level + 1);

	elements_n = 0;
	add_element(elements);
	elements.clear();
}

template<class element_type> element_type* octal_tree_node<element_type>::find_point(double x, double y, double z) {
	if(point_in_node(x,y,z)) {
		if(is_node_terminal) {
			for(int i = 0; i < elements_n; i++) 
				if (elements[i]->in_element(x,y,z))
					return elements[i];
		}
		else {
			for(int i = 0; i < 8; i++) {
				element_type* find = sub_nodes[i].find_point(x, y, z);
				if(find != NULL)
					return find;
			}

		}
	}
	return NULL;
}

template<class element_type> void octal_tree<element_type>::construct_tree(double x0, double x1, double y0, double y1, double z0, double z1, vector<element_type>& elements) {
	total_n = elements.size();
	split_constant = sqrt(total_n);
	int level_barier = log(total_n) / log(2.0) + 5; 
	root_node.init_node(x0, x1, y0, y1, z0, z1, split_constant, level_barier);
	root_node.set_level(0);
	vector<element_type*> el_v;
	for(int i = 0; i < total_n; i++) {
		el_v.push_back(&elements[i]);
	}
	root_node.add_element(el_v);
}

template<class element_type> element_type* octal_tree<element_type>::find_point(double x, double y, double z) {
	return root_node.find_point(x, y, z);
}