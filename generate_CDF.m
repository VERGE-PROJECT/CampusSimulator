function CDF_result = generate_CDF(vector_input)
%GENERATE_CDF This function generates the CDF of the vector and returns the
%result in two columns. The first column includes the values and the second
%column the probabilities.
    num_samples=size(vector_input,1);
    if num_samples>0
        histogram_size=min(num_samples,10000);
        CDF_result=zeros(histogram_size,2);
    
        [N,edges]=histcounts(vector_input,histogram_size);
    
        CDF_result(1,1)=0.5*(edges(1)+edges(2));
        CDF_result(1,2)=N(1)/num_samples;
        for i=2:histogram_size
            CDF_result(i,1)=0.5*(edges(i)+edges(i+1));
            CDF_result(i,2)=CDF_result(i-1,2)+N(i)/num_samples;
        end
    else
        CDF_result=[0,0];
    end
    
    %OLD FUNCTION
%     CDF_result=zeros(num_samples,2);
%     CDF_result(:,1)=sort(vector_input);
%     CDF_result(1,2)=1/num_samples;
%     for i=2:num_samples
%        CDF_result(i,2)=CDF_result(i-1,2)+1/num_samples;
%     end
end

