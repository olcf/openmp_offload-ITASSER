      subroutine metro(dE,id)
	use params
      common/ranzy/nozy
cc !$acc data copy(ranzy,nozy)

      id=1
      if(dE.gt.0)then
cc !$acc kernels
         if(aranzy(nozy).gt.weight12(dE))then
            id=3
         endif
cc !$acc end kernels
      endif
      return
cc !$acc end data
      end
