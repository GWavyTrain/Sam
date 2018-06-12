echo ''
echo Making binnable chunks for TrimmedAndCombinedData.dat
echo ''

Tot=7722072936
echo Total number of mergers: $Tot

ChunkSize=100000000
echo Chunk size: $ChunkSize

echo Total number of chunks: $(( $Tot / $ChunkSize + 1 )) +/- 1

chunk=1

nComments=13
Min=$(( $nComments + 1 ))
Max=$(( $Min + ( $ChunkSize - 1 ) ))

echo ''
echo chunk $chunk: $Min, $Max
#echo "sed -n '$(($Min)),$(($Max))p' TrimmedAndCombinedData.dat > ChunkedMergersUnique/mergers_chunk$chunk.dat"
sed -n "$(($Min)),$(($Max))p" TrimmedAndCombinedData.dat > ChunkedMergersUnique/mergers_chunk$chunk.dat

while [ $Max -lt $(( $Tot + $nComments )) ]
do
  chunk=$(( $chunk + 1 ))
  Min=$(( $Max + 1 ))
  Max=$(( $Min + ( $ChunkSize - 1 ) ))
  if [ $Max -gt $(( $Tot + $nComments )) ]
  then
    Max=$(( $Tot + $nComments ))
  fi
  echo chunk $chunk: $Min, $Max
#  echo "sed -n '$((Min)),$((Max))p' TrimmedAndCombinedData.dat > ChunkedMergersUnique/mergers_chunk$chunk.dat"
  sed -n "$((Min)),$((Max))p" TrimmedAndCombinedData.dat > ChunkedMergersUnique/mergers_chunk$chunk.dat
done

