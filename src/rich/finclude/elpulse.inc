
      INTEGER  nelcoo,typho,NELMAX
      INTEGER  ihitcc                   !   1900
      REAL  pulse
      PARAMETER (NELMAX=100000)
      INTEGER  kdophgen
c
cc1900      COMMON / ELPULSE / nelcoo(NELMAX,3),pulse(NELMAX,25),
      COMMON / ELPULSE / nelcoo(NELMAX,3),pulse(NELMAX,25),
     +                   typho(NELMAX),ihitcc(NELMAX),kdophgen
