#

set ii=0

while ($ii<100)
  @ ii++
  qsub <<FIN
  TestSpeed 5000 ~/results_TestSpeed.txt
FIN

end



