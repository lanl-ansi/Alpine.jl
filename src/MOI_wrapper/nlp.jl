"""
MOI NLPBlock 
"""
	
MOI.supports(::Optimizer, ::MOI.NLPBlock) = true
	
function MOI.set(model::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    model.nlp_data = nlp_data
    if model.nlp_data.has_objective 
        model.objective = nothing 
    end
end
