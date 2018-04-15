#CFLAGS	     = -Wall -pipe -O3 -march=prescott -fomit-frame-pointer\
#               -mfpmath=sse\
#               -ffast-math\
#               -I/usr/local/include
#LDFLAGS = -lm -L/usr/local/lib -lgmp
CFLAGS	     = -Wall -O3\
               -I/opt/local/include
LDFLAGS = -lm -L/opt/local/lib -lgmp

all:	two-body check-two-body-table


two_body_OBJ = two-body.o

two-body:	$(two_body_OBJ)
	$(CC) -o two-body $(two_body_OBJ) $(LDFLAGS)

check_two_body_table_OBJ = check-two-body-table.o

check-two-body-table:	$(check_two_body_table_OBJ)
	$(CC) -o check-two-body-table $(check_two_body_table_OBJ) $(LDFLAGS)

clean:
	$(RM) *.o *~ \
	two-body \
	check-two-body-table

clean-tables:
	$(RM) two-body.[XYZxyz]_[pqv]* two-body.TQ_[pqv]
