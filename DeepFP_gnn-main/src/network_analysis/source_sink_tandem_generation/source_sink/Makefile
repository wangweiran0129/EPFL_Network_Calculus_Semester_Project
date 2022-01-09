OUT=source_sink_pb2.py source_sink.descr

all: $(OUT)

%_pb2.py: %.proto
	protoc -I=$(dir $<) --python_out=$(dir $@) $<

%.descr: %.proto
	protoc -I=$(dir $<) --include_imports --descriptor_set_out=$@ $<

clean:
	rm -rf $(OUT)