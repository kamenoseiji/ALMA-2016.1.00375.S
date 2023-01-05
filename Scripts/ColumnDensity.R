# Column density using equivalent width
#-------- constants
k <- 1.380649e-16	# Boltzmann constant, erg K^-1
h <- 6.62607015e-27 # Planck constant, erg s
epsilon_0 <- 8.8541878128e-12 # F m^-1
h_k <- 4.799244662211e-11 # h / k (s K^-1)

ColumnFact <- function(mu, B, J, nu, T, EW){
	# mu : dipole moment (Debye)
	# B  : rotational frequency [MHz]
	# J  : lower angular momentum level
	# nu : frequency [MHz]
	# T  : excitation temperature [K]
	# EW : equivalent width = integ tau dv [km/s]
	hnu_kT <- h_k* nu* 1e6 / T	# h nu / k T
	hB_kT  <- h_k* B* 1e6 / T
	exp1 <- exp(hB_kT / 3.0)
	exp2 <- exp(hB_kT * J * (J+1))
	exp3 <- 1.0 / (exp(hnu_kT) - 1.0)
	Factor <- 3.0* k* T / (8.0* pi^3* (mu * 1.0e-18)^2 * B* 1.0e6 * (2* J + 1))
	return( Factor* exp1* exp2* exp3* EW*1e5 )
}

#-------- Sub-mm absorption
#CO21 <- ColumnFact(0.1101, 57636, 1, 230538, 230, c(17.6, 0.5))
# CO32 <- ColumnFact(0.1101, 57636, 2, 345790, 230, c(23.1, 0.5))
SO87_77 <- ColumnFact(1.535, 21523.556, 7, 214357, 344, c(0.9, 0.1))
SO55_44 <- ColumnFact(1.535, 21523.556, 4, 215221, 344, c(4.7, 0.1))
SO87_76 <- ColumnFact(1.535, 21523.556, 7, 340714, 344, c(8.1, 0.5))
SO88_77 <- ColumnFact(1.535, 21523.556, 7, 344311, 344, c(7.9, 0.5))
SO89_78 <- ColumnFact(1.535, 21523.556, 7, 346528, 344, c(11.7, 0.6))
S34O56_45 <- ColumnFact(1.535, 21102.73, 4, 215839, 344, c(0.9, 0.1))
CS76 <- ColumnFact(1.957, 24495.56, 6, 342883, 344, c(6.7, 0.3))
HC15N <- ColumnFact(2.9852, 43027.64, 3, 344200, 344, c(1.2, 0.6))
H13CN <- ColumnFact(2.9852, 43170.127, 3, 345340, 344, c(11.9, 0.5))
HCN43 <- ColumnFact(2.984, 44315.97, 3, 354505, 344, c(27.0,0.5))
HCNv2 <- ColumnFact(2.9420, 44422.42, 3, 356256, 344, c(13.0, 0.4))
HCO43 <- ColumnFact(3.888, 44594.4, 3, 356734, 344, c(15.9, 0.4))
HCOl1e <- ColumnFact(3.90, 44677.15, 3, 356549, 344, c(4.5, 0.4))
HCOl1f <- ColumnFact(3.90, 44677.15, 3, 358242, 344, c(5.3, 0.2))
#-------- mm absorption
SO22_11 <- ColumnFact(1.535, 21523.556, 1, 86083.95, 26, c(1.28, 0.07))
SO32_21 <- ColumnFact(1.535, 21523.556, 2, 99299.87, 26, c(2.33, 0.11))
SO33_22 <- ColumnFact(1.535, 21523.556, 2,129138.92, 26, c(3.02, 0.11))
HCN10 <- ColumnFact(2.984, 44315.97, 0, 88631.85, 26, c(4.95,0.07))
H13CN10 <- ColumnFact(2.984, 43170.127, 0, 86339.92, 26, c(0.64,0.10))
HCO10 <- ColumnFact(3.888, 44594.4, 0, 89188.53, 26, c(1.83, 0.13))
CS21 <- ColumnFact(1.957, 24495.56, 1, 97980.95, 26, c(0.47,0.13))

fcov <- 0.17
CO32corr  <- fcov * 44 /0.1508*0.1284 # correction factor for optically thick CO32
HCN43corr <- fcov * 44                # correction factor for optically thick HCN43
H13CNcorr <- fcov * 44 /0.1508*0.0594 # correction factor for optically thick H13CN
HCO43corr <- fcov * 44 /0.1508*0.0980 # correction factor for optically thick HCO+43

ColumnFact(0.1101, 57636, 2, 345790, 344, CO32corr* c(23.1, 0.5))
ColumnFact(2.984, 44315.97, 3, 354505, 344, HCN43corr* c(27.0,0.5))
ColumnFact(2.9852, 43170.127, 3, 345340, 344, H13CNcorr* c(11.9, 0.5))
ColumnFact(3.888, 44594.4, 3, 356734, 344, HCO43corr* c(15.9, 0.4))

#-------- average SO (mm) column density
mmSO <- c(SO22_11[1], SO32_21[1], SO33_22[1])
emmSO <- c(SO22_11[2], SO32_21[2], SO33_22[2])
weight <- 1/emmSO^2
mmSOstat <- c((mmSO %*% weight)/sum(weight), 1/sqrt(sum(weight)))
text_sd <- sprintf('N(SO, mm) : %.2e +- %.2e\n',  mmSOstat[1], mmSOstat[2] )
cat(text_sd)

#-------- average SO (mm) column density
submmSO <- c(SO89_78[1], SO88_77[1], SO87_76[1], SO55_44[1])
esubmmSO <- c(SO89_78[2], SO88_77[2], SO87_76[2], SO55_44[2])
weight <- 1/esubmmSO^2
submmSOstat <- c((submmSO %*% weight)/sum(weight), 1/sqrt(sum(weight)))
text_sd <- sprintf('N(SO, submm) : %.2e +- %.2e\n',  submmSOstat[1], submmSOstat[2] )
cat(text_sd)

#-------- Other sources from literatures
IRAS20551_HCO <- ColumnFact(3.888, 44594.4, 3, 356734, 36, c(16.9,0.2))     # IRAS 20551-4250 HCO+ J=4-3 (Imanishi+2017)
IRAS20551_SO  <- ColumnFact(1.535, 21523.556, 7, 344311, 36, c(0.64, 0.15)) # IRAS 20551-4250 SO 8_8 - 7_7   (Imanishi+2017)