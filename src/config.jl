abstract POD

type PODModel <: POD
	d :: JuMP.NLPEvaluator
	m :: JuMP.Model
	terms :: Array{Expr}
	t2x :: Dict  # Mapping of terms -> pod_x
	t2y :: Dict  # Mapping of terms -> pod_y
	PODModel() = new()
end


type PODConfig <: POD
	UBSolver :: String
	LBSolver :: String
	

end
