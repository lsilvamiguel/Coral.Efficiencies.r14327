rm -fr x
g++ -o x daq.cc -I ../xmlparse ../xmlparse/libexpat.a
x < daq.xml
