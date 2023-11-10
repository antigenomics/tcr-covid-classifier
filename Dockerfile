FROM avsastry/biopython:1.0
ENV port_num=5001
EXPOSE $port_num

COPY app docker_app/app
COPY .streamlit docker_app/.streamlit
COPY utils docker_app/utils
COPY figures/associations docker_app/figures/associations
COPY figures/*.csv docker_app/figures
COPY data/significant_clone_matrix_fisher_fmba_TRA_top_500k_wo_leaks.csv docker_app/data/
COPY data/significant_clone_matrix_fisher_fmba_TRB_top_500k_wo_leaks.csv docker_app/data/
COPY data/vdjdb.txt docker_app/data/
RUN pip install igraph
RUN pip install tqdm
RUN pip install streamlit
RUN pip install plotly
RUN pip install multipy
RUN pip install scipy
RUN pip install logomaker
WORKDIR docker_app
ENTRYPOINT streamlit run app/streamlit_app.py --server.port=$port_num
