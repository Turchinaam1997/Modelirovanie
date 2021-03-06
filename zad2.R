# �������� ������ � �������
dafr = read_csv2("city_trees2.csv")
# ������� �� �������� � �������
a = data.frame(dafr)
a
#������� ����� ������� ������ -- ��� ���� ����������
colnames(a)
# 1. ��� ���������� ����� ���������� ��� ������
# ��
# 2. ������� ������������� ����������
a$HR <- NULL
a$dbh_mm <- NULL
# 3. ������ ����������� �� ���� ����������
names(a)[names(a) == "dbh_m"] = "dbh"
# � ��������� ���������� ����������� �� �����������
# 4. ���� ���������� ������ �������� �����������
library(units) #������������� �������������� �����
#������������� ������ ����������� 
units(a$Ht) = as_units("m")
units(a$dbh) = as_units("m")
units(a$Clearance) = as_units("m")
units(a$Crown_Depth) = as_units("m")
units(a$Total_NSEW_Radial_Crown_Spread) = as_units("m")
units(a$Average_Radial_Crown_spread) = as_units("m")
units(a$Crown_Diameter) = as_units("m")
units(a$Predicted_CD_comb_f) = as_units("m")
units(a$Predicted_CD) = as_units("m")
# 5. ���� �����-�� ���������� �������� ������� ������ ����������, ��� ������ ���� ������ � ��������� � ���� ������ � �������� ����������
# ������� ����� ���������� � ����������� �� � ������
a = a %>% mutate(error = Predicted_CD_comb_f - Crown_Diameter)
a = a %>% mutate(error2 = Predicted_CD - Crown_Diameter)
# ������� ������ ��������
a$Diference <- NULL
a$Difference_comb_f <- NULL

# 6, 7, 8.  �������������� ���������� ������ ���� ���������, ��������� ���������� �� ����� ������ ���� ������, ���� �������������� ���������� �������� �� �����������
#library(forcats)
#library(sf)
#a$Data_Set1
#a = a %>%
#  mutate(Data_Set = as.factor("1", "0")) %>%
#  mutate(Data_Set = fct_recode(Data_Set, Norwich = "1", Peterborough = "0"))
#a$Data_Set
#---��� �� ��������. ��������� ������ ������
a$Data_Set[a$Data_Set == "0"] = "Peterborough"
a$Data_Set[a$Data_Set == "1"] = "Norwich"
a$Data_Set
#---� ��� ��� ��������. ������ ����� �� �����, ��������� �� ���� ��� �� :)

# 10. ������� ���� �� ������
# maple - Acer platanoides, 
# Oak - Quercus robur,
# Silver birch - Betula pendula, 
# Sycamore - Platanus occidentalis
a$Species[a$Species == "Oak"] = "Quercus robur"
a$Species[a$Species == "Norway_maple"] = "Acer platanoides"
a$Species[a$Species == "Silver_Birch"] = "Betula pendula"
a$Species[a$Species == "Sycamore"] = "Platanus occidentalis"

# 9. ������ ���� ������� ���������� ���������(lat,lon) � ���������� ������� ���������(� ������ ����� ���������) � � WGS84
library(stringr)
a$Grid_Reference
coord1 = str_replace_all(a$Grid_Reference, ' ', '')
coord_north = str_trunc(coord1, 12, "right", ellipsis = "") %>% str_trunc(5, "left", ellipsis = "")
coord_north
coord_e = str_trunc(coord1, 7, "right", ellipsis = "") %>% str_trunc(5, "left", ellipsis = "")
quadr = str_trunc(coord1, 2, "right", ellipsis = "")
table = data.frame(as.integer(coord_e), as.integer(coord_north), quadr)
names(table) = c("E", "N", "Quadr")
table = na.exclude(table)
#----------
table = table %>% mutate("Easting_BC" = case_when(
  quadr == "TF" ~ E + 600000, 
))
table = table %>% mutate("Northing_BC" = case_when(
  quadr == "TF" ~ N + 300000, 
))
table = table %>% mutate("Easting_BC" = case_when(
  quadr == "TG" ~ E + 700000, 
))
table = table %>% mutate("Northing_BC" = case_when(
  quadr == "TG" ~ N + 300000, 
))
table = table %>% mutate("Easting_BC" = case_when(
  quadr == "TL" ~ E + 600000, 
))
table = table %>% mutate("Northing_BC" = case_when(
  quadr == "TL" ~ N + 200000, 
))
table = na.exclude(table)
#-----------
#����������� ���������� ��� ��������
table_WGS = table %>% st_as_sf(coords = c("Easting_BC", "Northing_BC"), crs = 27700) %>% st_transform(4326) %>% st_coordinates() %>% as.data.frame()
write.csv2(a, file = "city_trees3.csv")