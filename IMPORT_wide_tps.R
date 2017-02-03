# A recipe to read such a file that contains XY columns, the others being the 
# covariates (or $fac) data.frame

# ID	Tank	Treatment	Status	Length	Weight	Week	 X1	 Y1	 X2	 Y2	 X3	 Y3	 X4	 Y4	 X5	 Y5	 X6	 Y6	 X7	 Y7	 X8	 Y8	 X9	 Y9	 X10	 Y10	 X11	 Y11	 X12	 Y12	 X13	 Y13	 X14	 Y14	 X15	 Y15
# file_1.jpg	D	Control	Initial	63	2.9	1	12.68022	7.44055	12.31154	7.8204	12.58526	7.72544	12.96511	7.65282	13.67453	7.27297	13.25558	8.1779	15.39502	7.08305	14.97607	8.33431	16.42284	7.13891	15.76369	8.16673	17.18812	7.35676	17.02054	7.98798	17.99251	7.32883	18.02602	7.91536	18.25505	7.65282
# file_2.jpg	D	Control	Initial	69	4	1	13.12517	9.13103	12.77904	9.44968	13.03726	9.44968	13.46579	9.38375	14.2844	8.90577	13.77895	9.90019	16.1139	8.71348	15.74031	10.09248	17.14128	8.76293	16.6633	9.9771	17.95989	9.10356	17.88846	9.83975	18.8609	9.03214	18.84991	9.7024	19.05869	9.38375
# ...
# ...

# Dependencies
library(dplyr)
library(tidyr)
library(Momocs)

PATH = "~/Desktop/data.txt"
# read the "wide" tps
tps <- PATH %>% 
  read.table(h=T, sep="\t") %>% 
  as.tbl()

# ids of X/Y columns
xy_cols <- tps %>% 
  colnames %>% 
  grep("^(X|Y)[[:digit:]]+$", .) 

# only the X/Y columns
df_xy <- tps %>% select(xy_cols)

# for each row, create a shape by filling a matrix
coo <- lapply(1:nrow(df_xy), function(i) matrix(df_xy[i, ], ncol=2, byrow=TRUE))

# all non- X/Y columns in tps = the $fac data.frame
fac <- tps %>% select(-xy_cols)

# eventually creates the Ldk object
final <- Ldk(coo, fac=fac)

# then you can
final %>% fgProcrustes() %>% stack

