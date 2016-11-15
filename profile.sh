SORT_FLAG='-s time'

python -m cProfile $SORT_FLAG main.py 1 test/ 5000 500 1000 0 0 1 > profileNew.out
