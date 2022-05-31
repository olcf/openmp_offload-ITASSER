      subroutine get_bond(axx,ayy,azz,al)
	use params
      implicit integer(i-z)
      common/ranzy/nozy
      athita=acos(1.-2.*aranzy(nozy)) !thita angle in random, [0,pi]
      aphi=2.*3.1415926*aranzy(nozy) !phi angle in random, [0,2pi]
      axx=al*sin(athita)*cos(aphi)
      ayy=al*sin(athita)*sin(aphi)
      azz=al*cos(athita)
      return
      end
