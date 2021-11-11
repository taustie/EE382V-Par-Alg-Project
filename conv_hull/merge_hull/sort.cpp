vector<pair<int, int>> merge(int start, int end, vector<pair<int, int>> points)
{
    int mid = (start + end) / 2;
    int k = mid + 1;
    vector<pair<int, int>> temp;

    while(start <= mid && k <= end)
    {
        if(points[start].first < points[k].first)
        {
            temp.push_back(points[start]);
            start++;
        }
        else
        {
            temp.push_back(points[k]);
            k++;    
        }
    }

    if(start < mid)
    {
        #pragma omp parallel
        {
            #pragma omp for
            for()
        }
    }

    if(k < end)
    {

    }
    return temp;
}

vector<pair<int, int>> merge_sort_seq(vector<pair<int, int>> points, int start, int end)
{
    int mid = (start + end) / 2;


}

void merge_sort_par(vector<pair<int, int>> points, int start, int end, int threads)
{
    if(threads == 1)
    {
        merge_sort_seq(points, start, end);
        return;
    }

    int mid = (start + end) / 2;
    vector<pair<int, int>> left;
    vector<pair<int, int>> right;
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            merge_sort_par(points, start, mid, threads / 2);
        }
        #pragma omp section
        {
            merge_sort_par(points, mid+1, end, threads - threads/2);
        }
    }

    vector<pair<int, int>> merged =  merge(start, end, points);
    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0, j = start; j <= end; i++, j++)
        {
            points[j] = merged[i];
        }
    }
}