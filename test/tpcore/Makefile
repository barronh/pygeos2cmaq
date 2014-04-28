test: geos2cmaq.bcon.20060801.nc geos2cmaq.icon.20060801.nc test.bcon test.icon

test.bcon: geos2cmaq.bcon.20060801.nc geos2cmaq.bcon.20060801.check.nc
	diff -s $^

test.icon: geos2cmaq.icon.20060801.nc geos2cmaq.icon.20060801.check.nc
	diff -s $^

geos2cmaq.bcon.20060801.nc: config.bcon.yaml BC.TEST.20060801 tracerinfo.dat diaginfo.dat bc_profile_CB05.dat METBDY3D_Benchmark_TP
	python -m pygeos2cmaq $<

geos2cmaq.icon.20060801.nc: config.icon.yaml BC.TEST.20060801 tracerinfo.dat diaginfo.dat bc_profile_CB05.dat METCRO3D_Benchmark_TP
	python -m pygeos2cmaq $<

pygeos2cmaq_test.zip:
	curl -LO curl -LO ftp://data.as.essie.ufl.edu/pub/geos2cmaq/pygeos2cmaq_test.zip

BC.TEST.20060801:
	make pygeos2cmaq_test.zip
	unzip pygeos2cmaq_test.zip

geos2cmaq.icon.20060801.check.nc geos2cmaq.bcon.20060801.check.nc METBDY3D_Benchmark_TP METCRO3D_Benchmark_TP: BC.TEST.20060801

clean:
	rm -rf geos2cmaq.icon.20060801.nc geos2cmaq.bcon.20060801.nc *.png *.png.nc

distclean:
	make clean
	rm -rf geos2cmaq.icon.20060801.check.nc geos2cmaq.bcon.20060801.check.nc METBDY3D_Benchmark_TP METCRO3D_Benchmark_TP BC.TEST.20060801 pygeos2cmaq_test.zip