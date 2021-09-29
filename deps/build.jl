using Libdl

nautyver = "nauty27rc4"
nautydir = "$(@__DIR__)/$nautyver/"

run(`tar xfz $(@__DIR__)/$nautyver.tar.gz -C $(@__DIR__)`)
cd(()->run(`./configure --enable-tls`), nautydir)
cd(()->run(`make nauty.a CCOBJ='${CC} -fPIC -c ${CFLAGS} -o $@'`), nautydir)

nautyobj = nautydir .* ["traces.o", "gtools.o", "naututil.o", "nautinv.o", "gutil1.o", "gutil2.o", "gtnauty.o", "naugroup.o", "nautycliquer.o", "nauty.o", "nautil.o", "nausparse.o", "naugraph.o", "schreier.o", "naurng.o"]

run(`gcc -shared -fPIC -I $nautydir -o $(@__DIR__)/nauty.$dlext $(@__DIR__)/extra.c $nautyobj`)
