      subroutine diff(n,pp,dt)
      real pp(*)
      y1=pp(1)
      pp(1)=0.0
      do jk=2,n
          y2=pp(jk)
          pp(jk)=(pp(jk)-y1)/dt
          y1=y2
      enddo
      return
      end    
