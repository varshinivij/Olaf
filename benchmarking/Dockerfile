FROM python:3.11
WORKDIR /Olaf 
COPY . /Olaf
RUN pip install -r benchmarking/requirements.txt
RUN chmod +x benchmarking/create_benchmark_env.sh && ./benchmarking/create_benchmark_env.sh
CMD ["bash", "benchmarking/run_interactive.sh"]
EXPOSE 5000
