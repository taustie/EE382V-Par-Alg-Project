#include "quick_hull.hpp"
#include "parallel_funcs.hpp"

Point * global_point;

Point * get_empty(){
	Point * tmp = new Point;
	tmp->x = global_point->x;
	tmp->y = global_point->y;
	return tmp;
}

Dist_Info neutral_distance(){
	Dist_Info maxDistResult;
	maxDistResult.index = 0;
	maxDistResult.line_dist = 0.0;
	maxDistResult.p1_dist = 0.0;
	return maxDistResult;
}

Dist_Info Dist_Max_Compare(Dist_Info a, Dist_Info b){
	if(a.line_dist < b.line_dist){
		return b;
	}
	else if(a.line_dist > b.line_dist){
		return a;
	}
	else if(a.line_dist == b.line_dist){
		if(a.p1_dist < b.p1_dist){
			return a;
		}
		else{
			return b;
		}
	}

	return neutral_distance(); // should never reach here
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


void get_max_dist_parallel(std::vector<Point *> &left_points, Point** max, Point *p1, Point *p2){
	Line line = get_line(*p1, *p2);

	Dist_Info maxDistResult = neutral_distance();
	long long int index;

	#pragma omp declare reduction \
	(maxDist:Dist_Info:omp_out=Dist_Max_Compare(omp_out, omp_in)) \
	initializer(omp_priv = neutral_distance())

	#pragma omp parallel for reduction(maxDist:maxDistResult)
	for(index = 0; index < left_points.size(); index++){
		double distance = get_distance(line, *left_points.at(index));
		if(maxDistResult.line_dist < distance){
			maxDistResult.index = index;
			maxDistResult.line_dist = distance;
			maxDistResult.p1_dist = get_distance_point(*left_points.at(index), *p1);
		}
		else if(maxDistResult.line_dist == distance){
			double distance_p1 = get_distance_point(*left_points.at(index), *p1);
			if(distance_p1 < maxDistResult.p1_dist){
				maxDistResult.index = index;
				maxDistResult.p1_dist = distance_p1;
			}
		}
	}

	*max = left_points.at(maxDistResult.index);
}
