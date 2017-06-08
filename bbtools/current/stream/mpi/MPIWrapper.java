package stream.mpi;

import mpi.*;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import java.nio.ByteBuffer;

import align2.Shared;
import align2.ListNum;
import stream.Read;

/**
 * Wraps MPI class functions for access by other programs.
 * All MPI calls should go through this class.
 * It should also set Shared.MPI fields such as MPI_RANK.
 * 
 * @author Jonathan Rood
 * @date Dec 9, 2014
 *
 */
public class MPIWrapper {

	public static void mpiInit(String[] args) {
		if(Shared.USE_MPI && Shared.USE_CRISMPI) {
			if(verbose){System.out.println("Running MPI Init.");}
			try {
				MPI.Init(args);
				Shared.MPI_RANK=MPI.COMM_WORLD.getRank();
				Shared.MPI_NUM_RANKS=MPI.COMM_WORLD.getSize();
			} catch (MPIException e) {
				e.printStackTrace();
			}
		}
	}

	public static void mpiFinalize() {
		if(Shared.USE_MPI && Shared.USE_CRISMPI) {
			if(verbose){System.out.println("Running MPI Finalize.");}
			try {
				MPI.Finalize();
			} catch (MPIException e) {
				e.printStackTrace();
			}
		}
	}

	public static void broadcast(ListNum<Read> ln){
		try {
			byte[] b=serialize(ln);
			int[] bLength={b.length};
			if(verbose){System.err.println("MPI:      Master broadcasting message size of " + b.length + ".");}
			MPI.COMM_WORLD.bcast(bLength,1,MPI.INT,0);
			if(verbose){System.err.println("MPI:      Master broadcasting reads, blocking.");}
			MPI.COMM_WORLD.bcast(b,b.length,MPI.BYTE,0);
		} catch (MPIException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void unicast(ListNum<Read> ln, final int toRank){
		try {
			byte[] b=serialize(ln);
			int[] bLength={b.length};
			if(verbose){System.err.println("MPI:      Master unicasting message size of " + b.length + " to rank " + toRank + ".");}
			MPI.COMM_WORLD.send(bLength,1,MPI.INT,toRank,40);
			if(verbose){System.err.println("MPI:      Master unicasting reads, blocking.");}
			MPI.COMM_WORLD.send(b,b.length,MPI.BYTE,toRank,50);
		} catch (MPIException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void iUnicast(ListNum<Read> ln, final int toRank){
		try {
			byte[] b=serialize(ln);
			int[] bLength={b.length};
			ByteBuffer bbLength=ByteBuffer.allocateDirect(4);
			bbLength.putInt(bLength[0]);
			bbLength.clear();
			if(verbose){System.err.println("MPI:      Master unicasting message size of " + b.length + " to rank " + toRank + ".");}
			MPI.COMM_WORLD.iSend(bbLength,4,MPI.BYTE,toRank,40);
			ByteBuffer bb=ByteBuffer.allocateDirect(b.length);
			bb.put(b);
			bb.clear();
			if(verbose){System.err.println("MPI:      Master unicasting reads, non-blocking.");}
			MPI.COMM_WORLD.iSend(bb,b.length,MPI.BYTE,toRank,50);
		} catch (MPIException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void iBroadcast(ListNum<Read> ln){
		try {
			byte[] b=serialize(ln);
			int[] bLength={b.length};
			ByteBuffer bbLength=ByteBuffer.allocateDirect(4);
			bbLength.putInt(bLength[0]);
			bbLength.clear();
			if(verbose){System.err.println("MPI:      Master broadcasting message size of " + b.length + ".");}
			MPI.COMM_WORLD.iBcast(bbLength,4,MPI.BYTE,0);
			ByteBuffer bb=ByteBuffer.allocateDirect(b.length);
			bb.put(b);
			bb.clear();
			if(verbose){System.err.println("MPI:      Master broadcasting reads, non-blocking.");}
			MPI.COMM_WORLD.iBcast(bb,b.length,MPI.BYTE,0);
		} catch (MPIException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void broadcastPaired(boolean b){
		try {
			boolean[] isPaired={b};
			if(verbose){System.err.println("MPI:      Master broadcasting pairing status of " + isPaired[0] + ".");}
			MPI.COMM_WORLD.bcast(isPaired,1,MPI.BOOLEAN,0);
		} catch (MPIException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void broadcastKeepall(boolean b){
		try {
			boolean[] keepAll={b};
			if(verbose){System.err.println("MPI:      Master broadcasting keepAll status of " + keepAll[0] + ".");}
			MPI.COMM_WORLD.bcast(keepAll,1,MPI.BOOLEAN,0);
		} catch (MPIException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static ListNum<Read> listen(){
		int rank=Shared.MPI_RANK;
		boolean keepAll=Shared.MPI_KEEP_ALL;
		if(keepAll){
			try {
				int[] bLength={0};
				MPI.COMM_WORLD.bcast(bLength,1,MPI.INT,0);
				byte[] b=new byte[bLength[0]];
				if(verbose){System.err.println("MPI:      Slave " + rank + " knows message will be of size " + bLength[0]);}
				MPI.COMM_WORLD.bcast(b,bLength[0],MPI.BYTE,0);
				if(verbose){System.err.println("MPI:      Slave " + rank + " received reads from blocking broadcast.");}
				ListNum<Read> ln=(ListNum<Read>) deserialize(b);
				return ln;
			} catch (MPIException e) {
				e.printStackTrace();
			} catch (IOException e){
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}else{
			try {
				int[] bLength={0};
				MPI.COMM_WORLD.recv(bLength,1,MPI.INT,0,40);
				byte[] b=new byte[bLength[0]];
				if(verbose){System.err.println("MPI:      Slave " + rank + " knows message will be of size " + bLength[0]);}
				MPI.COMM_WORLD.recv(b,bLength[0],MPI.BYTE,0,50);
				if(verbose){System.err.println("MPI:      Slave " + rank + " received reads from blocking unicast.");}
				ListNum<Read> ln=(ListNum<Read>) deserialize(b);
				return ln;
			} catch (MPIException e) {
				e.printStackTrace();
			} catch (IOException e){
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}

		}
		return null;
	}
	
	public static ListNum<Read> iListen(){
		int rank=Shared.MPI_RANK;
		boolean keepAll=Shared.MPI_KEEP_ALL;
		if(keepAll){
			try {
				int[] bLength={0};
				ByteBuffer bbLength=ByteBuffer.allocateDirect(4);
				Request req=MPI.COMM_WORLD.iBcast(bbLength,4,MPI.BYTE,0);
				req.waitFor();
				bbLength.clear();
				bLength[0]=bbLength.getInt();
				if(verbose){System.err.println("MPI:      Slave " + rank + " knows message will be of size " + bLength[0]);}
				byte[] b=new byte[bLength[0]];
				ByteBuffer bb=ByteBuffer.allocateDirect(bLength[0]);
				req=MPI.COMM_WORLD.iBcast(bb,bLength[0],MPI.BYTE,0);
				req.waitFor();
				bb.clear();
				bb.get(b);
				if(verbose){System.err.println("MPI:      Slave " + rank + " received reads from non-blocking broadcast.");}
				ListNum<Read> ln=(ListNum<Read>) deserialize(b);
				return ln;
			} catch (MPIException e) {
				e.printStackTrace();
			} catch (IOException e){
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}else{
			try {
				int[] bLength={0};
				ByteBuffer bbLength=ByteBuffer.allocateDirect(4);
				Request req=MPI.COMM_WORLD.iRecv(bbLength,4,MPI.BYTE,0,40);
				req.waitFor();
				bbLength.clear();
				bLength[0]=bbLength.getInt();
				byte[] b=new byte[bLength[0]];
				ByteBuffer bb=ByteBuffer.allocateDirect(bLength[0]);
				if(verbose){System.err.println("MPI:      Slave " + rank + " knows message will be of size " + bLength[0]);}
				req=MPI.COMM_WORLD.iRecv(bb,bLength[0],MPI.BYTE,0,50);
				req.waitFor();
				bb.clear();
				bb.get(b);
				if(verbose){System.err.println("MPI:      Slave " + rank + " received reads from non-blocking unicast.");}
				ListNum<Read> ln=(ListNum<Read>) deserialize(b);
				return ln;
			} catch (MPIException e) {
				e.printStackTrace();
			} catch (IOException e){
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return null;
	}

	public static boolean listenPaired(){
		int rank=Shared.MPI_RANK;
		try {
			boolean[] isPaired={false};
			if(verbose){System.err.println("MPI:      Slave " + rank + " listening for pairing status.");}
			MPI.COMM_WORLD.bcast(isPaired,1,MPI.BOOLEAN,0);
			if(verbose){System.err.println("MPI:      Slave " + rank + " pairing status is " + isPaired[0] + ".");}
			return isPaired[0];
		} catch (MPIException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}

	public static boolean listenKeepall(){
		int rank=Shared.MPI_RANK;
		try {
			boolean[] keepAll={false};
			if(verbose){System.err.println("MPI:      Slave " + rank + " listening for keepAll status.");}
			MPI.COMM_WORLD.bcast(keepAll,1,MPI.BOOLEAN,0);
			if(verbose){System.err.println("MPI:      Slave " + rank + " keepAll status is " + keepAll[0] + ".");}
			return keepAll[0];
		} catch (MPIException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}

	private static byte[] serialize(Object obj) throws IOException {
		ByteArrayOutputStream b = new ByteArrayOutputStream();
		ObjectOutputStream o = new ObjectOutputStream(b);
		o.writeObject(obj);
		return b.toByteArray();
	}

	private static Object deserialize(byte[] bytes) throws IOException, ClassNotFoundException {
		ByteArrayInputStream b = new ByteArrayInputStream(bytes);
		ObjectInputStream o = new ObjectInputStream(b);
		return o.readObject();
	}

	private static boolean verbose=true;
	private int rank=-1;
	
}
