CC = gfortran

pareto: pareto.f90
	$(CC) pareto.f90 -o pareto 

clean:
	rm -f pareto 

.PHONY: clean
