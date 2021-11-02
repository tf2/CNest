ViteRbi <- function(data=NULL, states=c(-1,0,1), normalstate=1, emiss = c(-1, 1, 0, 1, 1, 1), trans=c(0.99, 0.01, 0, 0.005, 0.99, 0.005, 0, 0.01, 0.99), ep=2, tp=3, active=T) {
	if(is.null(data)) { data(test); data = test; }
	jumpy = NULL; normalstates = vector();
	u = unique(data[,1])
	for(x in 1:length(u)) {
		d = data[data[,1]==u[x],]	
		res <- .C("ViteRbi"
				,"data" = as.double(d[,3])
				,"states" = as.double(rep(normalstate, length(d[,1])))
				,"emissions" = as.double(emiss)
				,"transitions" = as.double(trans)
				,"dN" = as.integer(length(d[,3]))
				,"sN" = as.integer(length(states))
				,"eN" = as.integer(ep)
				,"tN" = as.integer(tp)
				,"PACKAGE" = "ViteRbi")
		jumpy = rbind(jumpy, cbind(d[,1:3], res$states))
		
		if(active) {
			par(mfrow=c(2,1))
			plot(d[,3], xlab="Index", ylab="Value", main="Data")
			plot(res$states, xlab="Index", ylab="State", main="Estimated States")
			print("hit return to continue")
			scan("")
		}
	}
	
	invisible(jumpy)
}

extract_calls <- function(data) {
	all_cnv_calls = NULL
	u = unique(data[,1])
	for(z in 1:length(u)) {
		cnv_calls = NULL
		jumpy = data[data[,1]==u[z],]
		start_index = 1
		start_chr = u[z]
		start_start = jumpy[start_index,2]
		start_state = jumpy[start_index,4]
		for(x in 2:length(jumpy[,1])) {
			if(jumpy[x,4]!=start_state) {
				call = data.frame(start_chr, start_start, jumpy[x-1,2], start_state, mean(jumpy[start_index:(x-1), 3]), length(jumpy[start_index:(x-1), 3]), start_index, x-1)
				colnames(call) = c("chr", "start", "stop", "state", "mean_lr2", "number_probes", "start_index", "stop_index")
				cnv_calls = rbind(cnv_calls, call)
				start_chr = jumpy[x,1]
				start_start = jumpy[x,2]
				start_state = jumpy[x,4]
				start_index = x
			}
			if(x==length(jumpy[,1])) {
				call = data.frame(start_chr, start_start, jumpy[x,2], start_state, mean(jumpy[start_index:x, 3]), length(jumpy[start_index:x, 3]), start_index, x)
				colnames(call) = c("chr", "start", "stop", "state", "mean_lr2", "number_probes", "start_index", "stop_index")
				cnv_calls = rbind(cnv_calls, call)
			}
		}
		all_cnv_calls = rbind(all_cnv_calls, cnv_calls)
		print(z)
	}
return(all_cnv_calls)
}
