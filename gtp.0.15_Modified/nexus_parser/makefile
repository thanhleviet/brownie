libmy_structures.a *.o: *.c *.h
	gcc -enable-objc-gc -c *.c
	ar -r libmy_structures.a *.o
	ranlib libmy_structures.a
#libtool -static -o libmy_structures.a *.o
clean:
	rm -f *.o *.a
