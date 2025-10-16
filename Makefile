CC = cc

install : pixel.a raycast.a renderer
	cp renderer $(HOME)/bin

pixel.a : pixel.c
	cc -O3 pixel.c -c -o pixel.a

raycast.a : raycast.c
	cc -O3 raycast.c -c -o raycast.a

renderer : render.c
	cc -O3 pixel.a raycast.a render.c -o renderer -lm
