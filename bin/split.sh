awk '
    {
        for( i = 4; i <= NF; i=i+10 )
        {
        print(i)
        fn = i ".txt";
        print($1, $2, $3, $(i))>>fn;
        close( fn );
        }
    }
' canoes.reads_new.txt

# fn = i ".txt";
# printf( "%s %s %s\n", $1, $2, $(i) ) >>fn;
# close( fn );
##printf( "%s %s %s %s\n", $1, $2, $3, $(i) ) >>fn;
