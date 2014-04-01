# Steven Martell

.plotNaturalMortality <- function( repObj )
{
	# repObj is a list object that contains the model output (ie. the report file)
	with(repObj, {
		
		matplot(yr, (M_tot), type="l", xlab="Year", ylab="Natural Mortality")
		grid()
		
	})
	
}