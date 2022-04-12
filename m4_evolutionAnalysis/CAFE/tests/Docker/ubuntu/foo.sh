wallswitch=`$CC -Wall 2>&1`

if echo "$wallswitch" | grep -i "Error-Unknown" >/dev/null ; then
        echo PGI
else
        echo GCC
fi

