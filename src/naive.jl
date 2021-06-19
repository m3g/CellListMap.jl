#
# Function that uses the naive algorithm, for testing
#
function map_naive!(f,output,x,box)
  @unpack sides, cutoff_sq = box
  for i in 1:length(x)-1
    xᵢ = x[i]
    for j in i+1:length(x)
      xⱼ = wrapone(x[j],sides,xᵢ)
      d2 = distance_sq(xᵢ,xⱼ)
      if d2 <= cutoff_sq
        output = f(xᵢ,xⱼ,i,j,d2,output)
      end
    end
  end
  return output
end

#
# Function that uses the naive algorithm, for testing
#
function map_naive_two!(f,output,x,y,box)
  @unpack sides, cutoff_sq = box
  for i in 1:length(x)
    xᵢ = x[i]
    for j in 1:length(y)
      yⱼ = wrapone(y[j],sides,xᵢ)
      d2 = distance_sq(xᵢ,yⱼ)
      if d2 <= cutoff_sq
        output = f(xᵢ,yⱼ,i,j,d2,output)
      end
    end
  end
  return output
end


