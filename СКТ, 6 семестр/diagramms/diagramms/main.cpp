#include "gen_grid.h"
#include "portret_s.h"
#include "2d_func.h"
#include <conio.h>

const int omp_threads = 4;

void main()
{
	//сетка
	generate_2unreg_grid();

	//портрет
	SLAE_port_gen gen_2D;	
	gen_2D.gen();

	omp_set_num_threads(omp_threads);

	//считываем все данные
	normal_field_2D NF;
	NF.read();
	printf("data_true\n");
	NF.move_data_true();
	printf("enddatatrue\n");
	NF.move_inverse();
	printf("end\n");
	_getch();
}
