# OUT=messages_pb2.py messages.descr
# OUT=source_sink_pb2.py source_sink.descr
OUT=addressbook_pb2.py addressbook.descr

all: $(OUT)

%_pb2.py: %.proto
	protoc -I=$(dir $<) --python_out=$(dir $@) $<

%.descr: %.proto
	protoc -I=$(dir $<) --include_imports --descriptor_set_out=$@ $<

clean:
	rm -rf $(OUT)