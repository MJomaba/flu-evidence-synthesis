/* Copyright 2013 Edwin van Leeuwen.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#ifndef BSON_STREAM_H
#define BSON_STREAM_H
#include<map>
#include "bson/bson.h"

namespace mongo {

template<class T>
void operator>>( const mongo::BSONElement &bel, T &t ) {
	bel.Val( t );
}

void operator>>( const mongo::BSONElement &bel, size_t &t );
inline void operator>>( const mongo::BSONElement &bel, size_t &t ) {
	long long cpy = bel.number();
	//bel.Val( cpy );
	if (cpy>=0)
		t = (size_t) cpy;
	else
		throw MsgAssertionException(0, "Trying to convert negative number to size_t");
}

template<class T>
void operator>>( const mongo::BSONObj &bobj, T &t ) {
	// If we have a object we probably want to pipe the object as a whole into 
	// class T. To do that we first turn it into an element (containing an sub object)
	// and then pipe that to our class T. I could not find a simple way of turning an
	// object into an element, so we'll wrap it up and then when we unwrap it we get
	// an element back
	mongo::BSONObjBuilder bbuild;
	bbuild << "a" << bobj;
	bbuild.obj()["a"] >> t;
}

void operator>>( const mongo::BSONElement &bel, double &t );
inline void operator>>( const mongo::BSONElement &bel, double &t ) {
	t = bel.Number();
}

template<class T>
void operator>>( const mongo::BSONElement &bel, std::vector<T> &v ) {
	v.clear();
	auto barr = bel.Array();
	for ( auto & bson_el : barr ) {
		// This will only work if T has an empty constructor T()
		// I am not aware of a more general way of doing this though
		T el;
		bson_el >> el;
		v.push_back( el );
	}
}

template<class T>
void operator>>( const mongo::BSONElement &bel, std::list<T> &v ) {
	v.clear();
	auto barr = bel.Array();
	for ( auto & bson_el : barr ) {
		// This will only work if T has an empty constructor T()
		// I am not aware of a more general way of doing this though
		T el;
		bson_el >> el;
		v.push_back( el );
	}
}

template<class T>
void operator>>( const mongo::BSONElement &bel, std::set<T> &v ) {
	v.clear();
	auto barr = bel.Array();
	for ( auto & bson_el : barr ) {
		// This will only work if T has an empty constructor T()
		// I am not aware of a more general way of doing this though
		T el;
		bson_el >> el;
		v.insert( el );
	}
}

template<class K, class V>
void operator>>( const mongo::BSONElement &bel, std::pair<K,V> &p ) {
	auto barr = bel.Array();
	barr[0] >> p.first;
	barr[1] >> p.second;
}

template<class K, class V>
void operator>>( const mongo::BSONElement &bel, std::map<K,V> &map ) {
	map.clear();
	auto barr = bel.Array();
	for ( auto & bson_el : barr ) {
		std::pair<K,V> el;
		bson_el >> el;
		map.insert( el );
	}
}

	
template<class V>
void operator>>( const mongo::BSONObj &bobj, std::map<char *,V> &map ) {
	map.clear();
	for ( mongo::BSONObj::iterator i = bobj.begin(); i.more(); ) {
		mongo::BSONElement el = i.next();
		V value;
		el >> value;
		map[el.fieldName()] = value;
	}
}

template<class V>
void operator>>( const mongo::BSONObj &bobj, std::map<std::string,V> &map ) {
	map.clear();
	for ( mongo::BSONObj::iterator i = bobj.begin(); i.more(); ) {
		mongo::BSONElement el = i.next();
		V value;
		el >> value;
		map[el.fieldName()] = value;
	}
}

	class BSONEmitter;

	class BSONValueEmitter {
		public:
			BSONValueEmitter( BSONEmitter *pEmitter );

			template<class T>
				BSONEmitter &append( const T &t );

			BSONEmitter &append( const double &t );
			BSONEmitter &append( const long long &t );
			BSONEmitter &append( const size_t &t );
			BSONEmitter &append( const bool &t );
			BSONEmitter &append( const int &t );
			BSONEmitter &append( const std::string &t );
			BSONEmitter &append( const BSONArray &t );
			BSONEmitter &append( const BSONObj &t );
			BSONEmitter &append( const OID &t );

			void endField( const char *name ) {
				fieldName = name;
				builder.endField( name );
			}

			BSONEmitter *pEmitter;
			BSONObjBuilderValueStream builder;
		protected:
			const char *fieldName;
			
	};

	/**
	 * \brief Define an emitter for BSONObjects
	 *
	 * Behaviour is partly based on BSONObjBuilder, but different enough that
	 * we need a separate class
	 */
	class BSONEmitter {
		public:
			BSONEmitter() : builder( new BSONObjBuilder() ), v_emitter( this ) {}
			BSONEmitter( BSONObjBuilder *builder ) 
				: builder( builder ), v_emitter( this )
			{}

			BSONObj obj() {
				auto bobj = builder->obj();
				// This invalidates builder any way, so we can delete it
				delete builder;
				return bobj;
			}

			BSONValueEmitter &append( const std::string &name ) {
				v_emitter.builder.endField( name.c_str() );
				return v_emitter;
			}

			BSONValueEmitter &append( const char *name ) {
				v_emitter.endField( name );
				return v_emitter;
			}


			BSONObjBuilder *builder;
			BSONValueEmitter v_emitter;
	};

	inline BSONValueEmitter::BSONValueEmitter( BSONEmitter *pEmitter ) 
		: pEmitter( pEmitter ), builder( pEmitter->builder ) {
		}

	template<class T>
		BSONEmitter &BSONValueEmitter::append( const T &t ) {
			mongo::BSONEmitter b;
			b << t;
			return this->append( b.obj() );
		}

	inline BSONEmitter &BSONValueEmitter::append( const double &t ) {
		pEmitter->builder = &(builder << t);
		return (*pEmitter);
	}

	inline BSONEmitter &BSONValueEmitter::append( const long long &t ) {
		pEmitter->builder = &(builder << t);
		return (*pEmitter);
	}

	inline BSONEmitter &BSONValueEmitter::append( const size_t &t ) {
		// Casting to long long, which should be save enough
		long long cpy = (long long) t;
		pEmitter->builder = &(builder << cpy);
		return (*pEmitter);
	}

	inline BSONEmitter &BSONValueEmitter::append( const bool &t ) {
		pEmitter->builder = &(builder << t);
		return (*pEmitter);
	}

	inline BSONEmitter &BSONValueEmitter::append( const int &t ) {
		pEmitter->builder = &(builder << t);
		return (*pEmitter);
	}

	inline BSONEmitter &BSONValueEmitter::append( const std::string &t ) {
		pEmitter->builder = &(builder << t);
		return (*pEmitter);
	}

	inline BSONEmitter &BSONValueEmitter::append( const BSONArray &t ) {
		pEmitter->builder = &(builder << t);
		return (*pEmitter);
	}

	inline BSONEmitter &BSONValueEmitter::append( const BSONObj &t ) {
		pEmitter->builder = &(builder << t);
		return (*pEmitter);
	}

	inline BSONEmitter &BSONValueEmitter::append( const OID &t ) {
		pEmitter->builder->append( fieldName, t );
		return (*pEmitter);
	}


	class BSONArrayEmitter {
		public:
			BSONArrayEmitter() {}

			template<class T>
				BSONArrayEmitter &append( const T &t ) {
					mongo::BSONEmitter b;
					b << t;
					builder.append( b.obj() ); 
					return *this;
				}

			BSONArrayEmitter &append(	const double &t ) {
				builder.append( t );
				return *this;
			}

			BSONArrayEmitter &append(	const long long &t ) {
				builder.append( t );
				return *this;
			}

			BSONArrayEmitter &append(	const size_t &t ) {
				// Casting to long long, which should be save enough
				long long cpy = (long long) t;
				builder.append( cpy );
				return *this;
			}

			BSONArrayEmitter &append(	const bool &t ) {
				builder.append( t );
				return *this;
			}

			BSONArrayEmitter &append(	const int &t ) {
				builder.append( t );
				return *this;
			}

			BSONArrayEmitter &append(	const std::string &t ) {
				builder.append( t );
				return *this;
			}

			BSONArrayEmitter &append(	const BSONArray &t ) {
				builder.append( t );
				return *this;
			}

			BSONArrayEmitter &append(	const OID &t ) {
				builder.append( t );
				return *this;
			}

			BSONArray arr() {
				return builder.arr();
			}

			BSONArrayBuilder builder;
	};

	template<class V>
		BSONEmitter &operator<<( BSONEmitter &wrap, 
				const std::pair<const char *,V> &t ) {
			wrap << t.first << t.second;
			return wrap;
		}
	template<class V>
		BSONEmitter &operator<<( BSONEmitter &wrap, 
				const std::pair<const std::string,V> &t ) {
			wrap << t.first << t.second;
			return wrap;
		}

	template<class V>
		BSONEmitter &operator<<( BSONEmitter &wrap, 
				const std::map<char *,V> &t ) {
			for (auto & pair : t)
				wrap << pair;
			return wrap;
		}


	template<class V>
		BSONEmitter &operator<<( BSONEmitter &wrap, 
				const std::map<std::string,V> &t ) {
			for (auto & pair : t)
				wrap << pair;
			return wrap;
		}

		BSONEmitter &operator<<( BSONEmitter &wrap, 
				const OID &id );
		inline BSONEmitter &operator<<( BSONEmitter &wrap, const OID &id ) {
			wrap.builder->append( "_id", id );
			return wrap;
		}
	

template<class T>
mongo::BSONValueEmitter &operator<<( mongo::BSONEmitter &emitter, const T &t ) {
	return emitter.append( t );
} 

template<class T>
mongo::BSONEmitter &operator<<( mongo::BSONValueEmitter &emitter, const T &t ) {
	return emitter.append( t );
}

template<class T>
mongo::BSONArrayEmitter &operator<<( mongo::BSONArrayEmitter &barr,
		const T &t ) {
	return barr.append( t );
}


template<class T>
mongo::BSONEmitter &operator<<( mongo::BSONValueEmitter &bbuild, 
		const std::vector<T> &vt ) { 
	mongo::BSONArrayEmitter b;
	for ( const T &el : vt ) {
		b << el;
	}
	return bbuild.append( b.arr() );
}

template<class T>
mongo::BSONEmitter &operator<<( mongo::BSONValueEmitter &bbuild, 
		const std::set<T> &vt ) { 
	mongo::BSONArrayEmitter b;
	for ( const T &el : vt ) {
		b << el;
	}
	return bbuild.append( b.arr() );
}

template<class T>
mongo::BSONEmitter &operator<<( mongo::BSONValueEmitter &bbuild, 
		const std::list<T> &vt ) { 
	mongo::BSONArrayEmitter b;
	for ( const T &el : vt ) {
		b << el;
	}
	return bbuild.append( b.arr() );
}

template<class K, class V>
mongo::BSONEmitter &operator<<( mongo::BSONValueEmitter &bbuild, 
		const std::pair<K,V> &p ) { 
	mongo::BSONArrayEmitter b;
	b << p.first << p.second;
	return bbuild.append( b.arr() );
}

template<class K, class V>
mongo::BSONEmitter &operator<<( mongo::BSONValueEmitter &bbuild, 
		const std::map<K,V> &map ) { 
	mongo::BSONArrayEmitter b;
	for (auto &p : map) {
		mongo::BSONArrayEmitter b2;
		b2 << p.first << p.second;
		b.append( b2.arr() );
	}
	return bbuild.append( b.arr() );
}

template<class T>
mongo::BSONArrayEmitter &operator<<( mongo::BSONArrayEmitter &bbuild, 
		const std::vector<T> &vt ) { 
	mongo::BSONArrayEmitter b;
	for ( const T &el : vt ) {
		b << el;
	}
	return bbuild.append( b.arr() );
}

template<class T>
mongo::BSONArrayEmitter &operator<<( mongo::BSONArrayEmitter &bbuild, 
		const std::set<T> &vt ) { 
	mongo::BSONArrayEmitter b;
	for ( const T &el : vt ) {
		b << el;
	}
	return bbuild.append( b.arr() );
}

template<class T>
mongo::BSONArrayEmitter &operator<<( mongo::BSONArrayEmitter &bbuild, 
		const std::list<T> &vt ) { 
	mongo::BSONArrayEmitter b;
	for ( const T &el : vt ) {
		b << el;
	}
	return bbuild.append( b.arr() );
}

};

#endif
