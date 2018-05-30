echo ''
echo Making binnable chunks for mergers.dat
echo ''

Tot=7722072936
echo Total number of mergers: $Tot

ChunkSize=100000000
echo Chunk size: $ChunkSize

echo Total number of chunks: $(( Tot / ChunkSize + 1 ))

chunk=1

Min=2
Max=$(( $ChunkSize + $Min ))

echo ''
echo chunk $chunk: $Min, $Max
sed -n "$((Min)),$((Max))p" mergers.dat > ChunkedMergers/mergers_chunk$chunk.dat

while [ $Max -lt $Tot ]
do
  chunk=$(( $chunk + 1 ))
  Min=$(( $Max + 1 ))
  Max=$(( $Min + $ChunkSize ))
  if [ $Max -gt $Tot ]
  then
    Max=$Tot
  fi
  echo chunk $chunk: $Min, $Max
  sed -n "$((Min)),$((Max))p" mergers.dat > ChunkedMergers/mergers_chunk$chunk.dat
done

