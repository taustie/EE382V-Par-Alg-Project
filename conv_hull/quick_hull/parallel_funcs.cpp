#include "quick_hull.hpp"

Point * global_point;

Point * get_empty(){
	Point * tmp = new Point;
	tmp->x = global_point->x;
	tmp->y = global_point->y;
	return tmp;
}

Point * myMax(Point * a, Point * b){
    if(a->y > b->y){
    	return a;
    }
    else if(a->y < b->y){
    	return b;
    }
    else if(a->y == b->y){
    	if(a->x > b->x){
      		return a;
       	}
       	else{
       		return b;
       	}
    }

    return NULL; // this should never happen
}

Point * myMin(Point * a, Point * b){
    if(a->y > b->y){
    	return b;
    }
    else if(a->y < b->y){
    	return a;
    }
    else if(a->y == b->y){
    	if(a->x > b->x){
      		return b;
       	}
       	else{
       		return a;
       	}
    }

    return NULL; // this should never happen
}

void get_max_min_points_parallel(std::vector<Point *> &input_points, Point** max, Point** min){

	global_point = input_points.at(0);
	Point *tmp_min = input_points.at(0);
	Point *tmp_max = input_points.at(0);

    #pragma omp declare reduction \
        (maxPoint:Point *:omp_out=myMax(omp_out, omp_in)) \
        initializer(omp_priv = get_empty())

    #pragma omp declare reduction \
        (minPoint:Point *:omp_out=myMin(omp_out, omp_in)) \
        initializer(omp_priv = get_empty())

    #pragma omp parallel for reduction(maxPoint : tmp_max) reduction(minPoint:tmp_min)
    for(int i = 0; i < input_points.size(); i++){
        if(input_points.at(i)->y > tmp_max->y){
        	tmp_max = input_points.at(i);
        }
        else if(input_points.at(i)->y == tmp_max->y){
        	if(input_points.at(i)->x > tmp_max->x){
        		tmp_max = input_points.at(i);
        	}
        }

        if(input_points.at(i)->y < tmp_min->y){
        	tmp_min = input_points.at(i);
        }
        else if(input_points.at(i)->y == tmp_min->y){
        	if(input_points.at(i)->x < tmp_min->x){
        		tmp_min = input_points.at(i);
        	}
        }
    }

    *max = tmp_max;
    *min = tmp_min;

    std::cout << "max x = " << (*max)->x << std::endl;
    std::cout << "max y = " << (*max)->y << std::endl;
    std::cout << "min x = " << (*min)->x << std::endl;
    std::cout << "min y = " << (*min)->y << std::endl;
}
