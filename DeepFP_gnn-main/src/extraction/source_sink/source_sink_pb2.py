# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: source_sink.proto
"""Generated protocol buffer code."""
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor.FileDescriptor(
  name='source_sink.proto',
  package='netcal',
  syntax='proto3',
  serialized_options=None,
  create_key=_descriptor._internal_create_key,
  serialized_pb=b'\n\x11source_sink.proto\x12\x06netcal\"Q\n\x07Network\x12\n\n\x02id\x18\x01 \x01(\x05\x12\x1e\n\x06server\x18\x02 \x03(\x0b\x32\x0e.netcal.Server\x12\x1a\n\x04\x66low\x18\x03 \x03(\x0b\x32\x0c.netcal.Flow\"3\n\x06Server\x12\n\n\x02id\x18\x01 \x01(\x05\x12\x0c\n\x04rate\x18\x02 \x01(\x01\x12\x0f\n\x07latency\x18\x03 \x01(\x01\"\xc3\x01\n\x04\x46low\x12\n\n\x02id\x18\x01 \x01(\x05\x12\x0c\n\x04rate\x18\x02 \x01(\x01\x12\r\n\x05\x62urst\x18\x03 \x01(\x01\x12\x0c\n\x04path\x18\x04 \x03(\x05\x12\x1c\n\x04pmoo\x18\x05 \x01(\x0b\x32\x0e.netcal.Result\x12 \n\x06pmoofp\x18\x06 \x01(\x0b\x32\x10.netcal.FPResult\x12\x1f\n\x07\x64\x65\x62orah\x18\x07 \x01(\x0b\x32\x0e.netcal.Result\x12#\n\tdeborahfp\x18\x08 \x01(\x0b\x32\x10.netcal.FPResult\"\x1d\n\x06Result\x12\x13\n\x0b\x64\x65lay_bound\x18\x01 \x01(\x01\"T\n\x08\x46PResult\x12\x13\n\x0b\x64\x65lay_bound\x18\x01 \x01(\x01\x12\x33\n\x14\x65xplored_combination\x18\x02 \x03(\x0b\x32\x15.netcal.FPCombination\"\xa8\x01\n\rFPCombination\x12\x13\n\x0b\x64\x65lay_bound\x18\x01 \x01(\x01\x12H\n\x12\x66lows_prolongation\x18\x02 \x03(\x0b\x32,.netcal.FPCombination.FlowsProlongationEntry\x1a\x38\n\x16\x46lowsProlongationEntry\x12\x0b\n\x03key\x18\x01 \x01(\x05\x12\r\n\x05value\x18\x02 \x01(\x05:\x02\x38\x01\x62\x06proto3'
)




_NETWORK = _descriptor.Descriptor(
  name='Network',
  full_name='netcal.Network',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='id', full_name='netcal.Network.id', index=0,
      number=1, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='server', full_name='netcal.Network.server', index=1,
      number=2, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='flow', full_name='netcal.Network.flow', index=2,
      number=3, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=29,
  serialized_end=110,
)


_SERVER = _descriptor.Descriptor(
  name='Server',
  full_name='netcal.Server',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='id', full_name='netcal.Server.id', index=0,
      number=1, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='rate', full_name='netcal.Server.rate', index=1,
      number=2, type=1, cpp_type=5, label=1,
      has_default_value=False, default_value=float(0),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='latency', full_name='netcal.Server.latency', index=2,
      number=3, type=1, cpp_type=5, label=1,
      has_default_value=False, default_value=float(0),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=112,
  serialized_end=163,
)


_FLOW = _descriptor.Descriptor(
  name='Flow',
  full_name='netcal.Flow',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='id', full_name='netcal.Flow.id', index=0,
      number=1, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='rate', full_name='netcal.Flow.rate', index=1,
      number=2, type=1, cpp_type=5, label=1,
      has_default_value=False, default_value=float(0),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='burst', full_name='netcal.Flow.burst', index=2,
      number=3, type=1, cpp_type=5, label=1,
      has_default_value=False, default_value=float(0),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='path', full_name='netcal.Flow.path', index=3,
      number=4, type=5, cpp_type=1, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='pmoo', full_name='netcal.Flow.pmoo', index=4,
      number=5, type=11, cpp_type=10, label=1,
      has_default_value=False, default_value=None,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='pmoofp', full_name='netcal.Flow.pmoofp', index=5,
      number=6, type=11, cpp_type=10, label=1,
      has_default_value=False, default_value=None,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='deborah', full_name='netcal.Flow.deborah', index=6,
      number=7, type=11, cpp_type=10, label=1,
      has_default_value=False, default_value=None,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='deborahfp', full_name='netcal.Flow.deborahfp', index=7,
      number=8, type=11, cpp_type=10, label=1,
      has_default_value=False, default_value=None,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=166,
  serialized_end=361,
)


_RESULT = _descriptor.Descriptor(
  name='Result',
  full_name='netcal.Result',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='delay_bound', full_name='netcal.Result.delay_bound', index=0,
      number=1, type=1, cpp_type=5, label=1,
      has_default_value=False, default_value=float(0),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=363,
  serialized_end=392,
)


_FPRESULT = _descriptor.Descriptor(
  name='FPResult',
  full_name='netcal.FPResult',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='delay_bound', full_name='netcal.FPResult.delay_bound', index=0,
      number=1, type=1, cpp_type=5, label=1,
      has_default_value=False, default_value=float(0),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='explored_combination', full_name='netcal.FPResult.explored_combination', index=1,
      number=2, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=394,
  serialized_end=478,
)


_FPCOMBINATION_FLOWSPROLONGATIONENTRY = _descriptor.Descriptor(
  name='FlowsProlongationEntry',
  full_name='netcal.FPCombination.FlowsProlongationEntry',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='key', full_name='netcal.FPCombination.FlowsProlongationEntry.key', index=0,
      number=1, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='value', full_name='netcal.FPCombination.FlowsProlongationEntry.value', index=1,
      number=2, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=b'8\001',
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=593,
  serialized_end=649,
)

_FPCOMBINATION = _descriptor.Descriptor(
  name='FPCombination',
  full_name='netcal.FPCombination',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='delay_bound', full_name='netcal.FPCombination.delay_bound', index=0,
      number=1, type=1, cpp_type=5, label=1,
      has_default_value=False, default_value=float(0),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='flows_prolongation', full_name='netcal.FPCombination.flows_prolongation', index=1,
      number=2, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[_FPCOMBINATION_FLOWSPROLONGATIONENTRY, ],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=481,
  serialized_end=649,
)

_NETWORK.fields_by_name['server'].message_type = _SERVER
_NETWORK.fields_by_name['flow'].message_type = _FLOW
_FLOW.fields_by_name['pmoo'].message_type = _RESULT
_FLOW.fields_by_name['pmoofp'].message_type = _FPRESULT
_FLOW.fields_by_name['deborah'].message_type = _RESULT
_FLOW.fields_by_name['deborahfp'].message_type = _FPRESULT
_FPRESULT.fields_by_name['explored_combination'].message_type = _FPCOMBINATION
_FPCOMBINATION_FLOWSPROLONGATIONENTRY.containing_type = _FPCOMBINATION
_FPCOMBINATION.fields_by_name['flows_prolongation'].message_type = _FPCOMBINATION_FLOWSPROLONGATIONENTRY
DESCRIPTOR.message_types_by_name['Network'] = _NETWORK
DESCRIPTOR.message_types_by_name['Server'] = _SERVER
DESCRIPTOR.message_types_by_name['Flow'] = _FLOW
DESCRIPTOR.message_types_by_name['Result'] = _RESULT
DESCRIPTOR.message_types_by_name['FPResult'] = _FPRESULT
DESCRIPTOR.message_types_by_name['FPCombination'] = _FPCOMBINATION
_sym_db.RegisterFileDescriptor(DESCRIPTOR)

Network = _reflection.GeneratedProtocolMessageType('Network', (_message.Message,), {
  'DESCRIPTOR' : _NETWORK,
  '__module__' : 'source_sink_pb2'
  # @@protoc_insertion_point(class_scope:netcal.Network)
  })
_sym_db.RegisterMessage(Network)

Server = _reflection.GeneratedProtocolMessageType('Server', (_message.Message,), {
  'DESCRIPTOR' : _SERVER,
  '__module__' : 'source_sink_pb2'
  # @@protoc_insertion_point(class_scope:netcal.Server)
  })
_sym_db.RegisterMessage(Server)

Flow = _reflection.GeneratedProtocolMessageType('Flow', (_message.Message,), {
  'DESCRIPTOR' : _FLOW,
  '__module__' : 'source_sink_pb2'
  # @@protoc_insertion_point(class_scope:netcal.Flow)
  })
_sym_db.RegisterMessage(Flow)

Result = _reflection.GeneratedProtocolMessageType('Result', (_message.Message,), {
  'DESCRIPTOR' : _RESULT,
  '__module__' : 'source_sink_pb2'
  # @@protoc_insertion_point(class_scope:netcal.Result)
  })
_sym_db.RegisterMessage(Result)

FPResult = _reflection.GeneratedProtocolMessageType('FPResult', (_message.Message,), {
  'DESCRIPTOR' : _FPRESULT,
  '__module__' : 'source_sink_pb2'
  # @@protoc_insertion_point(class_scope:netcal.FPResult)
  })
_sym_db.RegisterMessage(FPResult)

FPCombination = _reflection.GeneratedProtocolMessageType('FPCombination', (_message.Message,), {

  'FlowsProlongationEntry' : _reflection.GeneratedProtocolMessageType('FlowsProlongationEntry', (_message.Message,), {
    'DESCRIPTOR' : _FPCOMBINATION_FLOWSPROLONGATIONENTRY,
    '__module__' : 'source_sink_pb2'
    # @@protoc_insertion_point(class_scope:netcal.FPCombination.FlowsProlongationEntry)
    })
  ,
  'DESCRIPTOR' : _FPCOMBINATION,
  '__module__' : 'source_sink_pb2'
  # @@protoc_insertion_point(class_scope:netcal.FPCombination)
  })
_sym_db.RegisterMessage(FPCombination)
_sym_db.RegisterMessage(FPCombination.FlowsProlongationEntry)


_FPCOMBINATION_FLOWSPROLONGATIONENTRY._options = None
# @@protoc_insertion_point(module_scope)
